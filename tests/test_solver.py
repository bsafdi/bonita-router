"""
Verification of the current-triangle solver and the field construction.

    python -m tests.test_solver      (or: pytest tests/)
"""
import datetime as dt
import numpy as np

from bonita.geo import BayGrid, OBSTRUCTIONS, SHIPPING_ZONES, MARKS, \
    DEFAULT_START_H, ll_to_xy, xy_to_ll, track_crosses_obstruction, \
    _segments_intersect
from bonita.polars import Polar
from bonita.router import Router, OFFSETS, route_course
from bonita.currents import build_current_field
from bonita.wind import build_wind_field


class _ZeroCurrent:
    def __init__(self, g): self.g = g
    def at(self, t): return np.zeros(self.g.X.shape), np.zeros(self.g.X.shape)


class _ConstWind:
    def __init__(self, g, kt, frm_deg):
        self.g, self.kt, self.frm = g, kt, np.deg2rad(frm_deg)
    def at(self, t):
        return np.full(self.g.X.shape, self.kt), np.full(self.g.X.shape, self.frm)


class _ConstCurrent:
    def __init__(self, g, u, v): self.g, self.u, self.v = g, u, v
    def at(self, t):
        return np.full(self.g.X.shape, self.u), np.full(self.g.X.shape, self.v)


def test_zero_current_matches_polar():
    """With no current, SOG on every grid heading must equal the polar exactly."""
    g = BayGrid(res_m=250.0)
    P = Polar()
    R = Router(g, _ConstWind(g, 14.0, 270.0), _ZeroCurrent(g), P)
    r, c = g.nearest_water_cell(-122.44, 37.83)
    _times, sog, _wd = R.edge_times(r, c, 12.0)
    for i, (dc, dr) in enumerate(OFFSETS):
        cog = np.degrees(np.arctan2(dc, dr)) % 360
        expect = float(P(cog - 270.0, 14.0))
        if expect <= R.min_sog:
            assert not np.isfinite(sog[i]) or sog[i] <= R.min_sog
        else:
            assert abs(sog[i] - expect) < 0.02, (cog, sog[i], expect)


def test_fair_current_adds_linearly():
    """Running dead downwind with a fair current: SOG = boatspeed + set."""
    g = BayGrid(res_m=250.0)
    P = Polar()
    R = Router(g, _ConstWind(g, 14.0, 270.0), _ConstCurrent(g, 2.0, 0.0), P)
    r, c = g.nearest_water_cell(-122.44, 37.83)
    _t, sog, _w = R.edge_times(r, c, 12.0)
    i = OFFSETS.index((1, 0))                       # COG 090
    assert abs(sog[i] - (float(P(180.0, 14.0)) + 2.0)) < 0.03


def test_impossible_course_is_infinite():
    """Dead upwind into 2 kt of foul current is not sailable."""
    g = BayGrid(res_m=250.0)
    R = Router(g, _ConstWind(g, 14.0, 270.0), _ConstCurrent(g, 2.0, 0.0), Polar())
    r, c = g.nearest_water_cell(-122.44, 37.83)
    times, _s, _w = R.edge_times(r, c, 12.0)
    assert not np.isfinite(times[OFFSETS.index((-1, 0))])


def test_fields_are_finite_over_the_race_window():
    g = BayGrid(res_m=250.0)
    d = dt.date(2026, 9, 12)
    cf = build_current_field(g, d)
    wf = build_wind_field(g, d)
    for t in np.arange(9.0, 21.0, 0.5):
        u, v = cf.at(t)
        s, dr = wf.at(t)
        for arr in (u, v, s, dr):
            assert np.all(np.isfinite(arr))
        assert np.nanmax(np.hypot(u, v)) < 6.0      # no runaway eddies


def test_polar_lands_near_the_rating():
    assert 120 < Polar().implied_phrf() < 185


# --------------------------------------------------------------------------
# SI 10 obstructions
# --------------------------------------------------------------------------

def _side_points(grid, pts, off_m):
    """
    For an obstruction polyline (metres), find a pair of navigable cells that
    straddle it: walk along the line and try perpendicular offsets of
    +-off_m at each sample.  Returns [(cellA, cellB), ...] for every sample
    where both sides are water.
    """
    pairs = []
    P = np.asarray(pts, float)
    for (ax, ay), (bx, by) in zip(P[:-1], P[1:]):
        L = np.hypot(bx - ax, by - ay)
        nx_, ny_ = -(by - ay) / L, (bx - ax) / L
        for f in np.linspace(0.05, 0.95, 19):
            mx, my = ax + f * (bx - ax), ay + f * (by - ay)
            cells = []
            for sgn in (+1, -1):
                lo, la = xy_to_ll(mx + sgn * off_m * nx_, my + sgn * off_m * ny_)
                rc = grid.nearest_water_cell(float(lo), float(la))
                cx, cy = grid.cell_xy(rc)
                # the nearest water cell must really be on that side
                if (cx - mx) * nx_ * sgn + (cy - my) * ny_ * sgn <= 0:
                    cells = None
                    break
                cells.append(rc)
            if cells:
                pairs.append(tuple(cells))
    return pairs


def test_every_obstruction_line_is_uncrossable():
    """
    For each SI 10 line: take cells straddling it and route across in a
    constant wind with no current.  Either there is no route, or the route
    goes round the end of the line -- it never crosses it.  Also every grid
    edge that straddles the line must be blocked outright.
    """
    g = BayGrid(res_m=250.0)
    R = Router(g, _ConstWind(g, 14.0, 250.0), _ZeroCurrent(g), Polar())
    checked = 0
    for (key, pts), ob in zip(g.obstruction_lines, OBSTRUCTIONS):
        assert ob["key"] == key
        # 1. direct edge test: every grid edge whose segment crosses the line
        #    is marked blocked, in both directions
        P = np.asarray(pts)
        for r in range(g.ny):
            for c in range(g.nx):
                for j, (dc, dr) in enumerate(OFFSETS):
                    rr, cc = r + dr, c + dc
                    if not (0 <= rr < g.ny and 0 <= cc < g.nx):
                        continue
                    ax, ay = g.X[r, c], g.Y[r, c]
                    bx, by = g.X[rr, cc], g.Y[rr, cc]
                    hit = any(_segments_intersect(ax, ay, bx, by, *p0, *p1)
                              for p0, p1 in zip(P[:-1], P[1:]))
                    if hit:
                        assert g.blocked[r, c, j], (key, r, c, j)
        # 2. routing test across the line
        pairs = _side_points(g, pts, off_m=g.obstruction_margin + 0.6 * g.res)
        for a, b in pairs[::3]:
            for src, dst in ((a, b), (b, a)):
                res = R.route(g.cell_center_ll(src), g.cell_center_ll(dst), 12.0)
                if not res.ok:
                    continue                      # sealed off entirely: fine
                xs, ys = ll_to_xy(res.lon, res.lat)
                assert not track_crosses_obstruction(
                    np.column_stack([xs, ys]), [(key, pts)]), key
                checked += 1
    assert checked > 0


def test_race_routes_clear_all_obstructions():
    """Both SI 12 courses, real fields, race-day start: no line crossed."""
    from bonita.currents import build_current_field
    from bonita.wind import build_wind_field
    g = BayGrid(res_m=250.0)
    d = dt.date(2026, 9, 12)
    R = Router(g, build_wind_field(g, d), build_current_field(g, d), Polar())
    for course in (1, 2):
        l1, l2 = route_course(R, DEFAULT_START_H, course=course)
        assert l1.ok and l2 is not None and l2.ok
        xs, ys = ll_to_xy(np.concatenate([l1.lon, l2.lon]),
                          np.concatenate([l1.lat, l2.lat]))
        assert not track_crosses_obstruction(np.column_stack([xs, ys]),
                                             g.obstruction_lines)
        # the Berkeley Pier daymark must be left to the west: every outbound
        # point south of the pier line must be west of the daymark
        ax, ay = ll_to_xy(*OBSTRUCTIONS[0]["line"][0])
        bx, by = ll_to_xy(*OBSTRUCTIONS[0]["line"][-1])
        x1, y1 = ll_to_xy(l1.lon, l1.lat)
        # A->B runs ENE, so a point south of the line has a positive cross product
        south = ((x1 - ax) * (by - ay) - (y1 - ay) * (bx - ax)) > 0
        assert np.all(x1[south] <= ax + 1e-6)


def test_berkeley_pier_daymark_position():
    """The pier obstruction starts at the ENC light and is > 2 NM long."""
    ob = OBSTRUCTIONS[0]
    assert ob["si"] == "10.1.a"
    lo, la = ob["line"][0]
    assert abs(la - 37.847735) < 2e-5 and abs(lo + 122.360551) < 2e-5
    xs, ys = ll_to_xy(*zip(*ob["line"]))
    L = np.sum(np.hypot(np.diff(xs), np.diff(ys))) / 1852.0
    assert 2.2 < L < 2.6


def test_time_along_track_refuses_a_crossing_track():
    from bonita.router import time_along_track
    g = BayGrid(res_m=250.0)
    R = Router(g, _ConstWind(g, 14.0, 250.0), _ZeroCurrent(g), Polar())
    # straight through the Fort Point passage: from Crissy to Baker Beach
    t = time_along_track(R, [(-122.470, 37.8095), (-122.4785, 37.8110),
                             (-122.486, 37.803)], 12.0)
    assert not np.isfinite(t)


def test_lane_cost_is_soft_and_flagged():
    """
    Lane cost must (a) never change real elapsed time bookkeeping -- the
    reported t_end is real time, (b) flag crossings, and (c) with weights
    off, route and cost coincide.
    """
    g = BayGrid(res_m=250.0)
    W = _ConstWind(g, 14.0, 250.0)
    C = _ZeroCurrent(g)
    a, b = (-122.4900, 37.8100), (-122.4400, 37.8300)   # Baker Beach -> N of Alcatraz
    r_soft = Router(g, W, C, Polar(), lane_cost=True).route(a, b, 12.0)
    r_flag = Router(g, W, C, Polar(), lane_cost=False).route(a, b, 12.0)
    assert r_soft.ok and r_flag.ok
    assert r_flag.lane_penalty_h == 0.0
    assert r_soft.lane_penalty_h >= 0.0
    assert r_soft.t_end >= r_flag.t_end - 1e-9        # soft cost cannot beat the free optimum
    assert any(z["kind"] == "lane" or z["kind"] == "precautionary"
               for z in r_flag.zone_crossings)         # that water is all zoned
    for z in r_soft.zone_crossings:
        assert z["t_in"] <= z["t_out"]
        assert z["key"] in {q["key"] for q in SHIPPING_ZONES}


def test_js_port_agrees_with_python():
    """
    The browser model must reproduce the Python grid (water, closed cells,
    blocked edges) and both race routes to within a few seconds.  Skipped
    when node is not installed.
    """
    import json, os, shutil, subprocess
    from bonita.currents import build_current_field
    from bonita.wind import build_wind_field
    node = shutil.which("node")
    if node is None:
        print("  (node not found -- JS parity not checked)")
        return
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    subprocess.run([shutil.which("python3") or "python", os.path.join(here, "export_model_data.py")],
                   check=True, capture_output=True)
    js = json.loads(subprocess.run(
        [node, os.path.join(here, "tests", "js_parity.js"),
         os.path.join(here, "web", "model_data.json")],
        check=True, capture_output=True, text=True).stdout)
    g = BayGrid(res_m=250.0)
    assert js["nWater"] == int(g.water.sum())
    assert js["nObstructed"] == int(g.obstructed.sum())
    assert js["nBlockedCells"] == int(g.blocked.any(axis=2).sum())
    d = dt.date(2026, 9, 12)
    R = Router(g, build_wind_field(g, d), build_current_field(g, d), Polar())
    for course in (1, 2):
        l1, l2 = route_course(R, DEFAULT_START_H, course=course)
        j = js["courses"][str(course)]
        assert abs(j["tRound"] - l1.t_end) * 3600 < 5, (course, j, l1.t_end)
        assert abs(j["tFinish"] - l2.t_end) * 3600 < 5, (course, j, l2.t_end)
        assert not j["crosses"]
        assert j["zonesOut"] == [z["key"] for z in l1.zone_crossings]
        assert j["zonesBack"] == [z["key"] for z in l2.zone_crossings]


def test_start_is_north_of_the_pier_and_navigable():
    g = BayGrid(res_m=250.0)
    r, c = g.nearest_water_cell(*MARKS["start"])
    assert g.water[r, c]
    sx, sy = ll_to_xy(*MARKS["start"])
    ax, ay = ll_to_xy(*OBSTRUCTIONS[0]["line"][0])
    bx, by = ll_to_xy(*OBSTRUCTIONS[0]["line"][-1])
    assert (sx - ax) * (by - ay) - (sy - ay) * (bx - ax) < 0   # north of the pier line


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
