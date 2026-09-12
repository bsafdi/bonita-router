"""
Strategy sweeps: what actually changes the answer.

Four questions this answers:

  1. Which side of the Bay?  Each named lane is a corridor; the router sails
     each one *optimally within that lane*, so the comparison is honest -- it
     is the cost of committing to a side, not the cost of steering badly.
  2. How much does the start time matter?  You cannot choose it, but it tells
     you how sharp the tidal gate is and so how expensive a bad start is.
  3. How much does the breeze scenario matter -- peak strength, and how early
     it fills?
  4. How much does the (constructed) polar matter?
"""

from __future__ import annotations

import numpy as np

from .geo import MARKS, ll_to_xy
from .router import Router, route_course
from .wind import SeaBreeze, build_wind_field
from .polars import Polar


# --------------------------------------------------------------------------
# Lanes
# --------------------------------------------------------------------------

# Corridors are keyed by course (geo.COURSES).  Course 1 rounds Point Bonita
# outside the Gate; course 2 rounds EASOM off Yellow Bluff, so its lanes are
# the three ways past Alcatraz.  Waypoints are lon/lat; the corridor is
# LANE_HALF_WIDTH_M either side of the polyline.  A route must stay inside
# the corridor AND clear of the SI 10 obstructions -- the Cityfront lane no
# longer includes the water inside Anita Rock or the Fort Point passage.
LANES = {
    1: {
        "out": {
            "Cityfront": [MARKS["start"], (-122.3900, 37.8330), (-122.4250, 37.8120),
                          (-122.4620, 37.8100), (-122.4780, 37.8160), MARKS["bonita"]],
            "mid-Bay (Alcatraz)": [MARKS["start"], (-122.3950, 37.8500),
                                   (-122.4240, 37.8330), (-122.4550, 37.8220),
                                   (-122.4780, 37.8200), MARKS["bonita"]],
            "Marin shore": [MARKS["start"], (-122.4050, 37.8600), (-122.4450, 37.8550),
                            (-122.4720, 37.8400), (-122.4790, 37.8280),
                            MARKS["bonita"]],
            "rhumb line": [MARKS["start"], MARKS["bonita"]],
        },
        "back": {
            "Marin shore": [MARKS["bonita"], (-122.4800, 37.8280), (-122.4760, 37.8390),
                            (-122.4600, 37.8520), (-122.4200, 37.8570),
                            MARKS["finish"]],
            "mid-Bay (Alcatraz)": [MARKS["bonita"], (-122.4780, 37.8200),
                                   (-122.4400, 37.8290), (-122.4050, 37.8480),
                                   MARKS["finish"]],
            "Cityfront": [MARKS["bonita"], (-122.4780, 37.8160), (-122.4400, 37.8110),
                          (-122.4000, 37.8250), (-122.3700, 37.8520), MARKS["finish"]],
            "rhumb line": [MARKS["bonita"], MARKS["finish"]],
        },
    },
    2: {
        "out": {
            "south of Alcatraz": [MARKS["start"], (-122.3900, 37.8330),
                                  (-122.4250, 37.8140), (-122.4450, 37.8180),
                                  (-122.4600, 37.8320), MARKS["easom"]],
            "north of Alcatraz": [MARKS["start"], (-122.3950, 37.8480),
                                  (-122.4250, 37.8360), (-122.4450, 37.8350),
                                  MARKS["easom"]],
            "Angel Island shore": [MARKS["start"], (-122.4050, 37.8600),
                                   (-122.4400, 37.8500), (-122.4550, 37.8430),
                                   MARKS["easom"]],
            "Raccoon Strait": [MARKS["start"], (-122.4000, 37.8780),
                               (-122.4350, 37.8760), (-122.4600, 37.8660),
                               MARKS["easom"]],
            "rhumb line": [MARKS["start"], MARKS["easom"]],
        },
        "back": {
            "Angel Island shore": [MARKS["easom"], (-122.4550, 37.8430),
                                   (-122.4400, 37.8500), (-122.4050, 37.8600),
                                   MARKS["finish"]],
            "north of Alcatraz": [MARKS["easom"], (-122.4450, 37.8350),
                                  (-122.4250, 37.8360), (-122.3950, 37.8480),
                                  MARKS["finish"]],
            "south of Alcatraz": [MARKS["easom"], (-122.4600, 37.8320),
                                  (-122.4450, 37.8180), (-122.4250, 37.8140),
                                  (-122.3900, 37.8330), MARKS["finish"]],
            "Raccoon Strait": [MARKS["easom"], (-122.4600, 37.8660),
                               (-122.4350, 37.8760), (-122.4000, 37.8780),
                               MARKS["finish"]],
            "rhumb line": [MARKS["easom"], MARKS["finish"]],
        },
    },
}

# Backward-compatible aliases (course 1).
OUTBOUND_LANES = LANES[1]["out"]
INBOUND_LANES = LANES[1]["back"]

LANE_HALF_WIDTH_M = 1500.0


def lane_mask(grid, waypoints_ll, half_width_m=LANE_HALF_WIDTH_M):
    """True where a cell is OUTSIDE the corridor (i.e. forbidden)."""
    dist = np.full(grid.X.shape, np.inf)
    for a, b in zip(waypoints_ll[:-1], waypoints_ll[1:]):
        ax, ay = ll_to_xy(*a)
        bx, by = ll_to_xy(*b)
        vx, vy = bx - ax, by - ay
        L2 = max(vx * vx + vy * vy, 1e-9)
        t = np.clip(((grid.X - ax) * vx + (grid.Y - ay) * vy) / L2, 0.0, 1.0)
        dist = np.minimum(dist, np.hypot(grid.X - (ax + t * vx),
                                         grid.Y - (ay + t * vy)))
    return dist > half_width_m


def compare_lanes(router, t_start, t_round=None, course=1,
                  half_width_m=LANE_HALF_WIDTH_M):
    """
    Route each lane optimally within its corridor.  Returns (outbound, inbound)
    lists of dicts sorted fastest first.  `t_round` fixes the departure time for
    the inbound comparison so the lanes are compared on equal terms.
    """
    from .geo import COURSES
    g = router.g
    mark = MARKS[COURSES[course]["mark"]]
    out = []
    for name, wps in LANES[course]["out"].items():
        res = router.route(MARKS["start"], mark, t_start,
                           forbid=lane_mask(g, wps, half_width_m))
        out.append(dict(name=name, ok=res.ok,
                        hours=res.elapsed_h if res.ok else np.inf,
                        arrive=res.t_end if res.ok else np.inf, route=res))
    out.sort(key=lambda r: r["hours"])

    t_round = t_round if t_round is not None else out[0]["arrive"]
    back = []
    for name, wps in LANES[course]["back"].items():
        res = router.route(mark, MARKS["finish"], t_round,
                           forbid=lane_mask(g, wps, half_width_m))
        back.append(dict(name=name, ok=res.ok,
                         hours=res.elapsed_h if res.ok else np.inf,
                         finish=res.t_end if res.ok else np.inf, route=res))
    back.sort(key=lambda r: r["hours"])
    return out, back


# --------------------------------------------------------------------------
# Sweeps
# --------------------------------------------------------------------------

def sweep_start_times(grid, wind_field, current_field, polar=None,
                      starts=np.arange(9.75, 11.51, 0.25), scenario="base",
                      course=1, **router_kw):
    polar = polar or Polar()
    rows = []
    # One Router for the whole sweep: the edge cache is keyed by (cell, time
    # bin), not by start time, so later solves reuse most of the earlier work.
    R = Router(grid, wind_field, current_field, polar, **router_kw)
    for t0 in starts:
        l1, l2 = route_course(R, float(t0), course=course)
        if l2 is None or not l2.ok:
            rows.append(dict(scenario=scenario, start=float(t0),
                             total_h=np.nan, t_round=np.nan, finish=np.nan))
            continue
        rows.append(dict(scenario=scenario, start=float(t0),
                         round_h=l1.elapsed_h, back_h=l2.elapsed_h,
                         total_h=l2.t_end - l1.t_start,
                         t_round=l1.t_end, finish=l2.t_end))
    return rows


def sweep_wind_scenarios(grid, current_field, date, polar=None, t_start=10.25,
                         scenarios=None, course=1, **router_kw):
    polar = polar or Polar()
    scenarios = scenarios or [
        ("soft (peak 13 kt)", dict(peak_kt=13.0)),
        ("expected (peak 19 kt)", dict(peak_kt=19.0)),
        ("windy (peak 25 kt)", dict(peak_kt=25.0)),
        ("late fill (onset 12:00)", dict(peak_kt=19.0, t_onset=12.0)),
        ("early fill (onset 09:00)", dict(peak_kt=19.0, t_onset=9.0)),
    ]
    rows = []
    for label, kw in scenarios:
        wf = build_wind_field(grid, date, breeze=SeaBreeze(**kw))
        R = Router(grid, wf, current_field, polar, **router_kw)
        l1, l2 = route_course(R, t_start, course=course)
        if l2 is None or not l2.ok:
            rows.append(dict(scenario=label, total_h=np.nan))
            continue
        rows.append(dict(scenario=label, start=t_start,
                         round_h=l1.elapsed_h, back_h=l2.elapsed_h,
                         total_h=l2.t_end - l1.t_start,
                         t_round=l1.t_end, finish=l2.t_end,
                         wind_field=wf, legs=(l1, l2)))
    return rows


def sweep_polar_scale(grid, wind_field, current_field, t_start=10.25,
                      scales=(0.92, 0.96, 1.0, 1.04, 1.08), course=1, **router_kw):
    rows = []
    shared_fields = {}
    for s in scales:
        p = Polar(scale=s)
        R = Router(grid, wind_field, current_field, p, **router_kw)
        R._field_cache = shared_fields      # wind/current do not depend on polar
        l1, l2 = route_course(R, t_start, course=course)
        rows.append(dict(scale=s, implied_phrf=p.implied_phrf(),
                         total_h=(l2.t_end - l1.t_start) if l2 and l2.ok else np.nan,
                         t_round=l1.t_end if l1.ok else np.nan))
    return rows
