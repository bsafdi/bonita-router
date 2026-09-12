"""
Time-optimal routing through a time-dependent wind and current field.

Method
------
Earliest-arrival Dijkstra on a graph whose nodes are (grid cell, direction of
arrival).  Carrying the arrival direction in the state is what lets the router
charge a realistic penalty for tacking and gybing, which in turn is what stops
a grid router from producing the classic physically meaningless sawtooth: with
the penalty in place a beat comes out as a handful of long tacks, as it should.

Each edge is priced by solving the current triangle properly.  The grid fixes
the course over ground; the solver searches boat headings for the one whose
water velocity, added to the tidal set, points along that COG, and takes the
resulting speed over ground.  If no heading can hold the course -- current
stronger than the boat -- the edge is impassable.

Because leaving later can never get you there earlier in this field, the
earliest-arrival labels satisfy the FIFO property and plain Dijkstra is exact
up to the time discretisation.
"""

from __future__ import annotations

import heapq
from dataclasses import dataclass, field

import numpy as np

from .polars import chop_factor

# 16 grid moves, roughly every 22.5 degrees.  (dc = east, dr = north)
OFFSETS = [(1, 0), (2, 1), (1, 1), (1, 2), (0, 1), (-1, 2), (-1, 1), (-2, 1),
           (-1, 0), (-2, -1), (-1, -1), (-1, -2), (0, -1), (1, -2), (1, -1),
           (2, -1)]
NDIR = len(OFFSETS)
START_DIR = NDIR          # virtual "no previous heading" state


@dataclass
class RouteResult:
    ok: bool
    t_start: float
    t_end: float
    path_rc: list = field(default_factory=list)
    path_t: list = field(default_factory=list)
    lon: np.ndarray = None
    lat: np.ndarray = None
    note: str = ""
    # Shipping-lane bookkeeping (soft cost, SI 9.1).  `cost_end` is the
    # penalised arrival "time" the search minimised; `t_end` is real time.
    cost_end: float = np.inf
    zone_crossings: list = field(default_factory=list)

    @property
    def elapsed_h(self):
        return self.t_end - self.t_start

    @property
    def lane_penalty_h(self):
        """Extra search cost charged for time spent inside shipping zones."""
        return (self.cost_end - self.t_end) if np.isfinite(self.cost_end) else 0.0


class Router:
    def __init__(self, grid, wind_field, current_field, polar,
                 time_bin_min=5.0, tack_penalty_s=14.0, gybe_penalty_s=9.0,
                 n_heading_samples=240, min_sog_kt=0.15, chop_k=0.045,
                 current_scale=1.0, max_edge_bins=3.0, lane_cost=True):
        self.g = grid
        # Soft shipping-lane cost (geo.SHIPPING_ZONES).  When on, the time
        # spent on an edge that ENTERS a zone cell is multiplied by that
        # zone's weight in the search objective only; real elapsed time is
        # tracked separately and is what gets reported.  lane_cost=False
        # keeps the crossing flags but charges nothing.
        self.lane_cost = lane_cost
        self.zone_w = (grid.zone_weight if lane_cost
                       else np.ones_like(grid.zone_weight))
        self.wf = wind_field
        self.cf = current_field
        self.polar = polar
        self.dt_bin = time_bin_min / 60.0
        self.tack_pen = tack_penalty_s / 3600.0
        self.gybe_pen = gybe_penalty_s / 3600.0
        self.min_sog = min_sog_kt
        self.chop_k = chop_k
        self.current_scale = current_scale
        # An edge is priced with the field frozen at its DEPARTURE bin.  Let an
        # edge run for many bins and the router starts exploiting stale water:
        # it will happily "sail" a 559 m cell at 0.18 kt for 101 minutes,
        # through a flood that has long since eased in reality.  Anything that
        # takes longer than a few bins is not a leg, it is waiting -- and this
        # formulation cannot represent waiting, so it is refused instead.
        self.max_edge_h = max_edge_bins * self.dt_bin

        # geometry of the 16 moves
        res = grid.res
        d = np.array([[dc * res, dr * res] for dc, dr in OFFSETS], float)
        self.leg_m = np.hypot(d[:, 0], d[:, 1])
        self.leg_nm = self.leg_m / 1852.0
        self.ghat = d / self.leg_m[:, None]                    # (16,2) unit E/N
        self.cog_deg = np.rad2deg(np.arctan2(d[:, 0], d[:, 1])) % 360.0

        self.theta = np.linspace(0, 2 * np.pi, n_heading_samples, endpoint=False)
        self.sin_t = np.sin(self.theta)
        self.cos_t = np.cos(self.theta)
        self.theta_deg = np.rad2deg(self.theta)

        self._field_cache = {}
        self._edge_cache = {}

    # -- fields -----------------------------------------------------------

    def fields(self, t):
        """Wind and current at the time-bin containing t (cached)."""
        b = int(np.floor(t / self.dt_bin))
        hit = self._field_cache.get(b)
        if hit is None:
            tb = (b + 0.5) * self.dt_bin
            tws, twd = self.wf.at(tb)
            cu, cv = self.cf.at(tb)
            hit = (np.asarray(tws), np.asarray(twd),
                   np.asarray(cu) * self.current_scale,
                   np.asarray(cv) * self.current_scale)
            self._field_cache[b] = hit
        return b, hit

    # -- the current triangle ---------------------------------------------

    def edge_times(self, r, c, t):
        """
        Hours to traverse each of the 16 outgoing edges from cell (r,c),
        departing at time t.  np.inf where the edge cannot be held.
        """
        b, (tws, twd, cu, cv) = self.fields(t)
        key = (r, c, b)
        hit = self._edge_cache.get(key)
        if hit is not None:
            return hit

        w = float(tws[r, c])
        wd = float(twd[r, c])
        cx, cy = float(cu[r, c]), float(cv[r, c])

        # boat velocity through the water for every candidate heading
        twa = np.rad2deg(self.theta) - np.rad2deg(wd)
        bs = self.polar(twa, w)                                 # (nθ,)
        # chop penalty depends on the course being sailed
        bs = bs * chop_factor(wd, w, cx, cy, self.theta_deg, k=self.chop_k)
        vx = bs * self.sin_t + cx                               # (nθ,) over ground
        vy = bs * self.cos_t + cy

        gx = self.ghat[:, 0][:, None]                           # (16,1)
        gy = self.ghat[:, 1][:, None]
        cross = vx[None, :] * gy - vy[None, :] * gx             # (16,nθ)
        dot = vx[None, :] * gx + vy[None, :] * gy

        cn = np.roll(cross, -1, axis=1)
        dn = np.roll(dot, -1, axis=1)
        sign_change = (cross <= 0) & (cn > 0) | (cross >= 0) & (cn < 0)
        denom = cross - cn
        f = np.where(np.abs(denom) > 1e-12, cross / np.where(
            np.abs(denom) > 1e-12, denom, 1.0), 0.0)
        f = np.clip(f, 0.0, 1.0)
        sog = dot * (1 - f) + dn * f                            # interpolated SOG
        sog = np.where(sign_change & (sog > self.min_sog), sog, -np.inf)
        best = sog.max(axis=1)                                  # (16,)

        times = np.where(best > 0, self.leg_nm / np.maximum(best, 1e-9), np.inf)
        times = np.where(times > self.max_edge_h, np.inf, times)
        self._edge_cache[key] = (times, best, wd)
        return self._edge_cache[key]

    # -- turn cost --------------------------------------------------------

    def turn_penalty(self, i_in, i_out, wd):
        if i_in == START_DIR or i_in == i_out:
            return 0.0
        a_in = ((self.cog_deg[i_in] - np.rad2deg(wd) + 180.0) % 360.0) - 180.0
        a_out = ((self.cog_deg[i_out] - np.rad2deg(wd) + 180.0) % 360.0) - 180.0
        if a_in * a_out < 0:                    # crossed the wind axis
            # Which axis it crossed is what tells a tack from a gybe: the sum
            # of the two angles from the wind is < 180 if you went through
            # head-to-wind, > 180 if you went through dead downwind.
            return self.tack_pen if abs(a_in) + abs(a_out) < 180 else self.gybe_pen
        return 0.0

    # -- the search -------------------------------------------------------

    def route(self, start_ll, end_ll, t_start, t_limit=None, forbid=None):
        """
        Earliest-arrival search from start_ll to end_ll departing at t_start.

        The heap key is the *cost*: real elapsed time plus the soft shipping-
        lane surcharge.  Real arrival time is carried alongside so the fields
        are always sampled at the true clock time and the reported finish is
        the true finish.  With lane_cost=False (or all weights 1.0) cost and
        time coincide and this is exact earliest-arrival Dijkstra.

        SI 10 obstructions are hard: edges in `grid.blocked` never relax, and
        cells inside the obstruction margin are not water.
        """
        g = self.g
        ny, nx = g.ny, g.nx
        r0, c0 = g.nearest_water_cell(*start_ll)
        r1, c1 = g.nearest_water_cell(*end_ll)
        water = g.water if forbid is None else (g.water & ~forbid)
        if not water[r0, c0] or not water[r1, c1]:
            return RouteResult(False, t_start, np.inf,
                               note="start or end cell is not navigable")
        blocked = g.blocked
        zone = g.zone
        zw = self.zone_w

        NS = NDIR + 1
        best = np.full((ny, nx, NS), np.inf)          # cost labels
        treal = np.full((ny, nx, NS), np.inf)         # real time at label
        prev = np.full((ny, nx, NS, 3), -1, dtype=np.int32)
        best[r0, c0, START_DIR] = t_start
        treal[r0, c0, START_DIR] = t_start
        heap = [(t_start, r0, c0, START_DIR)]
        goal = None
        t_cap = np.inf if t_limit is None else t_start + t_limit

        while heap:
            k, r, c, di = heapq.heappop(heap)
            if k > best[r, c, di] + 1e-12:
                continue
            t = treal[r, c, di]
            if (r, c) == (r1, c1):
                goal = (r, c, di, t, k)
                break
            if t > t_cap:
                continue
            times, _sog, wd = self.edge_times(r, c, t)
            brow = blocked[r, c]
            for j, (dc, dr) in enumerate(OFFSETS):
                if brow[j] or not np.isfinite(times[j]):
                    continue
                rr, cc = r + dr, c + dc
                if not (0 <= rr < ny and 0 <= cc < nx) or not water[rr, cc]:
                    continue
                dt = times[j] + self.turn_penalty(di, j, wd)
                nt = t + dt
                nk = k + dt * zw[zone[rr, cc]]
                if nk < best[rr, cc, j] - 1e-12:
                    best[rr, cc, j] = nk
                    treal[rr, cc, j] = nt
                    prev[rr, cc, j] = (r, c, di)
                    heapq.heappush(heap, (nk, rr, cc, j))

        if goal is None:
            return RouteResult(False, t_start, np.inf, note="no route found")

        r, c, di, t_end, k_end = goal
        path, times = [], []
        while True:
            path.append((r, c))
            times.append(float(treal[r, c, di]))
            pr, pc, pdi = prev[r, c, di]
            if pr < 0:
                break
            r, c, di = int(pr), int(pc), int(pdi)
        path.reverse()
        times.reverse()
        lon = np.array([g.LON[a, b] for a, b in path])
        lat = np.array([g.LAT[a, b] for a, b in path])
        res = RouteResult(True, t_start, float(t_end), path, times, lon, lat,
                          cost_end=float(k_end))
        res.zone_crossings = zone_crossings(g, path, times)
        return res


def zone_crossings(grid, path_rc, path_t):
    """
    Where and when a routed path is inside a shipping zone.  Returns a list
    of dicts {zone, key, kind, t_in, t_out} in path order; consecutive cells
    in the same zone are merged into one crossing.
    """
    from .geo import SHIPPING_ZONES
    out = []
    cur = None
    for (r, c), t in zip(path_rc, path_t):
        k = int(grid.zone[r, c])
        if cur is not None and k == cur["_k"]:
            cur["t_out"] = t
            continue
        if cur is not None:
            out.append(cur)
            cur = None
        if k >= 0:
            z = SHIPPING_ZONES[k]
            cur = dict(_k=k, zone=z["name"], key=z["key"], kind=z["kind"],
                       t_in=t, t_out=t)
    if cur is not None:
        out.append(cur)
    for o in out:
        o.pop("_k")
    return out


def route_course(router, t_start, marks=None, course=1, **kw):
    """
    Route the whole race: start -> the course's turning mark -> finish,
    chaining the legs so the second one departs when the first one arrives.
    `course` is a key of geo.COURSES (1 = Point Bonita, 2 = EASOM).
    """
    from .geo import MARKS, COURSES
    marks = marks or MARKS
    mark = COURSES[course]["mark"]
    leg1 = router.route(marks["start"], marks[mark], t_start, **kw)
    if not leg1.ok:
        return leg1, None
    leg2 = router.route(marks[mark], marks["finish"], leg1.t_end, **kw)
    return leg1, leg2


# --------------------------------------------------------------------------
# Forced-route comparator: price a track you choose, with the same physics
# --------------------------------------------------------------------------

def time_along_track(router, waypoints_ll, t_start, step_m=None):
    """
    Elapsed time to sail a prescribed track (a list of lon/lat waypoints).

    NOTE: this maximises velocity made good along the track direction without
    requiring the boat to hold the track, so on a beat it is an optimistic
    LOWER BOUND, not a realistic time.  For comparing strategies, use
    `sweep.compare_lanes`, which routes each option properly inside a corridor.
    """
    from .geo import ll_to_xy
    g = router.g
    step = step_m or g.res
    pts = []
    for (lo0, la0), (lo1, la1) in zip(waypoints_ll[:-1], waypoints_ll[1:]):
        x0, y0 = ll_to_xy(lo0, la0)
        x1, y1 = ll_to_xy(lo1, la1)
        n = max(int(np.hypot(x1 - x0, y1 - y0) / step), 1)
        for k in range(n):
            f = k / n
            pts.append((x0 + f * (x1 - x0), y0 + f * (y1 - y0)))
    pts.append(ll_to_xy(*waypoints_ll[-1]))

    # a prescribed track that crosses an SI 10 obstruction is not a track
    from .geo import track_crosses_obstruction
    if track_crosses_obstruction(pts, g.obstruction_lines):
        return np.inf

    t = t_start
    x0g, y0g = g.x[0], g.y[0]
    for (ax, ay), (bx, by) in zip(pts[:-1], pts[1:]):
        dx, dy = bx - ax, by - ay
        L = np.hypot(dx, dy)
        if L < 1e-6:
            continue
        gx, gy = dx / L, dy / L
        ci = int(np.clip(round((ax - x0g) / g.res), 0, g.nx - 1))
        ri = int(np.clip(round((ay - y0g) / g.res), 0, g.ny - 1))
        if not g.water[ri, ci]:
            return np.inf
        _b, (tws, twd, cu, cv) = router.fields(t)
        w, wd = float(tws[ri, ci]), float(twd[ri, ci])
        # NOTE: router.fields() has already applied current_scale.
        cx, cy = float(cu[ri, ci]), float(cv[ri, ci])

        twa = router.theta_deg - np.rad2deg(wd)
        bs = router.polar(twa, w) * chop_factor(
            wd, w, cx, cy, router.theta_deg, k=router.chop_k)
        vx = bs * router.sin_t + cx
        vy = bs * router.cos_t + cy
        # velocity made good along the track direction, maximised over heading
        vmg = vx * gx + vy * gy
        best = float(np.max(vmg))
        if best <= router.min_sog:
            return np.inf
        t += (L / 1852.0) / best
    return t - t_start
