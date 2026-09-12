#!/usr/bin/env python3
"""
BYC Big Windward-Leeward Race to Point Bonita -- optimal route calculator.

Examples
--------
    # everything, with the defaults (cached NOAA currents + parameterised breeze)
    python run_race.py

    # race week: pull the live forecast and the SFBOFS current model
    python run_race.py --wind openmeteo --currents sfbofs

    # a specific scenario
    python run_race.py --start 10:20 --peak-kt 22 --onset 11:00

    # just the sensitivity study
    python run_race.py --only sweeps
"""

from __future__ import annotations

import argparse
import datetime as dt
import os

import numpy as np

from bonita.geo import BayGrid, COURSES, MARKS, DEFAULT_START_H, TIME_LIMIT_H, \
    OBSTRUCTIONS, obstruction_lines_ll
from bonita.currents import build_current_field
from bonita.wind import build_wind_field, SeaBreeze
from bonita.polars import Polar
from bonita.router import Router, route_course
from bonita.plotting import plot_route, plot_sweep, leg_table, \
    track_length_nm, hhmm
from bonita.sweep import compare_lanes, sweep_start_times, \
    sweep_wind_scenarios, sweep_polar_scale

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def parse_hhmm(s):
    if ":" in s:
        h, m = s.split(":")
        return int(h) + int(m) / 60.0
    return float(s)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--date", default="2026-09-12")
    ap.add_argument("--start", default=None,
                    help="start-gun time, local (default: SI 3.1 warning 1100 "
                         "+ 5 min = 11:05)")
    ap.add_argument("--course", default="both", choices=["1", "2", "both"],
                    help="SI 12.0 course: 1 = Point Bonita, 2 = EASOM, or both")
    ap.add_argument("--no-obstructions", action="store_true",
                    help="ignore the SI 10 obstruction lines (for comparison only)")
    ap.add_argument("--no-lane-cost", action="store_true",
                    help="flag shipping-lane crossings but do not penalise them")
    ap.add_argument("--skip-unconstrained", action="store_true",
                    help="do not also solve the unconstrained route for comparison")
    ap.add_argument("--res", type=float, default=250.0,
                    help="grid resolution in metres (200 is finer/slower)")
    ap.add_argument("--currents", default="station",
                    choices=["station", "sfbofs"])
    ap.add_argument("--wind", default="param", choices=["param", "openmeteo"])
    ap.add_argument("--peak-kt", type=float, default=19.0,
                    help="afternoon sea-breeze peak, domain mean")
    ap.add_argument("--base-kt", type=float, default=8.0)
    ap.add_argument("--onset", default="10:00", help="hour the breeze fills")
    ap.add_argument("--wind-peak-time", default="15:30")
    ap.add_argument("--polar-scale", type=float, default=1.0)
    ap.add_argument("--current-scale", type=float, default=1.0,
                    help="multiply the whole current field (spring/neap fudge)")
    ap.add_argument("--avoid-shoal", action="store_true",
                    help="treat the Berkeley flats as unnavigable")
    ap.add_argument("--only", default="all",
                    choices=["all", "route", "sweeps"])
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    date = dt.date.fromisoformat(args.date)
    t0 = parse_hhmm(args.start) if args.start else DEFAULT_START_H
    courses = [1, 2] if args.course == "both" else [int(args.course)]

    print(f"== BYC Big Windward-Leeward -- {date} start {hhmm(t0)} "
          f"(time limit {hhmm(TIME_LIMIT_H)}) ==\n")

    grid = BayGrid(res_m=args.res, avoid_shoal=args.avoid_shoal,
                   obstructions=not args.no_obstructions)
    print(f"grid {grid.ny}x{grid.nx} @ {args.res:.0f} m, "
          f"{int(grid.water.sum())} navigable cells")
    if not args.no_obstructions:
        print("SI 10 obstructions (hard):")
        for ob, (_k, pts) in zip(OBSTRUCTIONS, obstruction_lines_ll()):
            lo, la = pts[0]
            print(f"  {ob['si']:7s} {ob['name'][:58]:58s} {la:.5f}N {-lo:.5f}W  "
                  f"[{ob['status'].split(' ')[0]}]")
        print(f"  {int(grid.obstructed.sum())} cells closed, "
              f"{int(grid.blocked.any(axis=2).sum())} cells with blocked edges")
    print()

    breeze = SeaBreeze(base_kt=args.base_kt, peak_kt=args.peak_kt,
                       t_onset=parse_hhmm(args.onset),
                       t_peak=parse_hhmm(args.wind_peak_time))
    cf = build_current_field(grid, date, source=args.currents)
    wf = build_wind_field(grid, date, source=args.wind, breeze=breeze)
    polar = Polar(scale=args.polar_scale)
    print(f"polar: {polar.name}, scale {polar.scale:.2f}, "
          f"implied PHRF ~{polar.implied_phrf():.0f}")
    base = cf if hasattr(cf, "n_measured") else getattr(cf, "fallback", None)
    if base is not None:
        from bonita.currents import REF_STATION
        ref = base.series.get(REF_STATION) or next(iter(base.series.values()))
        print(f"current field: {base.n_measured} NOAA stations + "
              f"{len(base.zones) - base.n_measured} derived zones")
        print("reference station (kt, + = flood): " + "  ".join(
            f"{hhmm(h)} {float(ref(h)):+.1f}" for h in range(10, 20)))
    print()

    R = Router(grid, wf, cf, polar, current_scale=args.current_scale,
               lane_cost=not args.no_lane_cost)

    # Unconstrained comparison: same physics, no obstruction lines, no lane
    # cost.  This is what the model recommended before SI 10 was applied.
    R_free = None
    if args.only in ("all", "route") and not args.skip_unconstrained \
            and not args.no_obstructions:
        g_free = BayGrid(res_m=args.res, avoid_shoal=args.avoid_shoal,
                         obstructions=False)
        cf_free = build_current_field(g_free, date, source=args.currents)
        wf_free = build_wind_field(g_free, date, source=args.wind, breeze=breeze)
        R_free = Router(g_free, wf_free, cf_free, polar,
                        current_scale=args.current_scale, lane_cost=False)

    def crossings_str(leg):
        if not leg.zone_crossings:
            return "none"
        return "; ".join(f"{z['key']} {hhmm(z['t_in'])}-{hhmm(z['t_out'])}"
                         for z in leg.zone_crossings)

    if args.only in ("all", "route"):
        for course in courses:
            cinfo = COURSES[course]
            mark = cinfo["mark"]
            print(f"================ COURSE {course}: start -> {cinfo['name']} "
                  f"-> finish (~{cinfo['approx_nm']} NM) ================")
            l1, l2 = route_course(R, t0, course=course)
            if not (l1.ok and l2 and l2.ok):
                print("no route found:", l1.note if not l1.ok else l2.note)
                print()
                continue
            total = l2.t_end - l1.t_start
            makes_it = l2.t_end <= TIME_LIMIT_H
            print(f"  start          {hhmm(l1.t_start)}")
            print(f"  round {mark:8s} {hhmm(l1.t_end)}   "
                  f"(leg 1 {l1.elapsed_h:.2f} h, {track_length_nm(l1):.1f} nm sailed)")
            print(f"  finish         {hhmm(l2.t_end)}   "
                  f"(leg 2 {l2.elapsed_h:.2f} h, {track_length_nm(l2):.1f} nm sailed)")
            print(f"  ELAPSED        {total:.2f} h  ({int(total)}h{int(round((total%1)*60)):02d})")
            margin = (TIME_LIMIT_H - l2.t_end) * 60
            print(f"  TIME LIMIT     {hhmm(TIME_LIMIT_H)}  -> "
                  f"{'OK, ' + f'{margin:.0f} min in hand' if makes_it else 'DNF -- misses the limit by ' + f'{-margin:.0f} min'}")
            print(f"  shipping zones out : {crossings_str(l1)}")
            print(f"  shipping zones back: {crossings_str(l2)}")
            pen = (l1.lane_penalty_h + l2.lane_penalty_h) * 60
            if pen > 0.05:
                print(f"  (lane surcharge carried in the search: {pen:.1f} min; "
                      f"not part of the elapsed time above)")

            if R_free is not None:
                f1, f2 = route_course(R_free, t0, course=course)
                if f1.ok and f2 and f2.ok:
                    ftot = f2.t_end - f1.t_start
                    print(f"  unconstrained  round {hhmm(f1.t_end)}  finish "
                          f"{hhmm(f2.t_end)}  elapsed {ftot:.2f} h  -> SI 10 "
                          f"obstructions + lane cost add {(total - ftot)*60:+.1f} min")
                    from bonita.geo import track_crosses_obstruction, ll_to_xy
                    xs, ys = ll_to_xy(np.concatenate([f1.lon, f2.lon]),
                                      np.concatenate([f1.lat, f2.lat]))
                    crosses = track_crosses_obstruction(
                        np.column_stack([xs, ys]), grid.obstruction_lines)
                    print(f"  unconstrained route crosses an SI 10 line: "
                          f"{'YES -- that route would mean retiring' if crosses else 'no'}")
            print()

            print("  waypoint log (every 20 min)")
            for tag, leg in (("out ", l1), ("back", l2)):
                for r in leg_table(R, leg, 20):
                    print(f"    {tag} {r['time']}  {r['lat']:.4f} {r['lon']:.4f}"
                          f"   TWS {r['tws']:4.1f} kt from {r['twd']:3.0f}"
                          f"   current {r['cur_kt']:4.2f} kt -> {r['cur_set']:3.0f}")
            print()

            out_rows, back_rows = compare_lanes(R, t0, t_round=l1.t_end,
                                                course=course)
            print("  outbound: cost of committing to one side "
                  "(each lane sailed optimally within a 1.5 km corridor)")
            for r in out_rows:
                if not r["ok"]:
                    print(f"    {r['name']:22s}  no route inside that corridor")
                    continue
                d = (r["hours"] - l1.elapsed_h) * 60
                print(f"    {r['name']:22s} {r['hours']:5.2f} h   "
                      f"{d:+6.1f} min vs unconstrained optimum")
            print("  inbound (all leaving the mark at the same time)")
            for r in back_rows:
                if not r["ok"]:
                    print(f"    {r['name']:22s}  no route inside that corridor")
                    continue
                d = (r["hours"] - l2.elapsed_h) * 60
                print(f"    {r['name']:22s} {r['hours']:5.2f} h   "
                      f"{d:+6.1f} min vs unconstrained optimum")
            print()

            tag = f"c{course}"
            plot_route(grid, cf, wf, [l1, l2],
                       t_snapshots=[t0 + 0.5, t0 + 1.75, l1.t_end,
                                    min(l2.t_end, l1.t_end + 1.0)],
                       title=f"Course {course} ({cinfo['name']}) over tidal "
                             f"current  |  {date} start {hhmm(t0)}  |  "
                             f"elapsed {total:.2f} h",
                       path=os.path.join(OUT, f"route_current_{tag}.png"),
                       show="current")
            plot_route(grid, cf, wf, [l1, l2],
                       t_snapshots=[t0 + 1.0, l1.t_end + 0.5],
                       title=f"Course {course} over wind field",
                       path=os.path.join(OUT, f"route_wind_{tag}.png"), show="wind")
            print(f"  charts -> {OUT}/route_current_{tag}.png, "
                  f"{OUT}/route_wind_{tag}.png\n")

    if args.only in ("all", "sweeps"):
        for course in courses:
            print(f"---------------- SENSITIVITY (course {course}) ----------------")
            rows = sweep_start_times(grid, wf, cf, polar,
                                     starts=np.arange(10.5, 12.01, 0.25),
                                     course=course)
            print("  start time")
            for r in rows:
                print(f"    start {hhmm(r['start'])}  round {hhmm(r['t_round'])}"
                      f"  finish {hhmm(r['finish'])}  elapsed {r['total_h']:.2f} h")
            plot_sweep(rows, path=os.path.join(OUT, f"sweep_start_c{course}.png"))

            print("\n  breeze scenario (start fixed)")
            for r in sweep_wind_scenarios(grid, cf, date, polar, t_start=t0,
                                          course=course):
                if np.isnan(r["total_h"]):
                    print(f"    {r['scenario']:26s} no route")
                    continue
                print(f"    {r['scenario']:26s} out {r['round_h']:4.2f} h  "
                      f"back {r['back_h']:4.2f} h  elapsed {r['total_h']:5.2f} h  "
                      f"finish {hhmm(r['finish'])}")

            print("\n  polar scale (how much the constructed polar matters)")
            for r in sweep_polar_scale(grid, wf, cf, t_start=t0, course=course):
                print(f"    scale {r['scale']:.2f}  implied PHRF "
                      f"{r['implied_phrf']:5.0f}   elapsed {r['total_h']:5.2f} h")
            print(f"\n  chart -> {OUT}/sweep_start_c{course}.png\n")


if __name__ == "__main__":
    main()
