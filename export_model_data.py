#!/usr/bin/env python3
"""
Export the Python model's tables to web/model_data.json.

`bonita/` is the source of truth for geometry, current stations, derived zones,
wind structure and the polar; this script is the only thing that copies them
into the browser app.  Run it (or build_web.py, which calls it) after editing
any table.
"""
import datetime as dt
import json
import pathlib

from bonita import currents, geo, polars, sweep, wind

HERE = pathlib.Path(__file__).parent


def main():
    data = {
        "meta": {"date": "2026-09-12",
                 "race": "BYC Big Windward-Leeward Regatta (PB / EASOM)",
                 "tz": "America/Los_Angeles",
                 "built": dt.datetime.utcnow().isoformat() + "Z"},
        "domain": geo.DOMAIN,
        "proj": {"lat0": geo.LAT0, "lon0": geo.LON0,
                 "mPerDegLat": geo.M_PER_DEG_LAT,
                 "mPerDegLon": geo.M_PER_DEG_LON},
        "marks": {k: list(v) for k, v in geo.MARKS.items()},
        "land": [{"name": n, "height": geo.LAND_HEIGHT_M.get(n, 50.0),
                  "poly": [[round(a, 5), round(b, 5)] for a, b in p]}
                 for n, p in geo.LAND_POLYGONS.items()],
        # SI 12.0 courses and SI 3.1 / 6.2 times
        "race": {"startH": geo.DEFAULT_START_H, "warningH": geo.FIRST_WARNING_H,
                 "timeLimitH": geo.TIME_LIMIT_H},
        "courses": {str(k): v for k, v in geo.COURSES.items()},
        "yraMarks": {"DOC": list(geo.YRA_DOC), "FOC": list(geo.YRA_FOC)},
        # SI 10 obstruction lines, already finished at the closest shore and
        # sealed 200 m inland (geo.obstruction_lines_ll) -- hard barriers.
        "obstructionMarginM": geo.OBSTRUCTION_MARGIN_M,
        "obstructions": [{"key": ob["key"], "si": ob["si"], "name": ob["name"],
                          "status": ob["status"],
                          "line": [[round(a, 6), round(b, 6)] for a, b in pts]}
                         for ob, (_k, pts) in zip(geo.OBSTRUCTIONS,
                                                  geo.obstruction_lines_ll())],
        # SI 9.1 shipping lanes / precautionary areas -- soft cost + flag.
        "shippingZones": [{"key": z["key"], "kind": z["kind"], "name": z["name"],
                           "weight": z["weight"],
                           "poly": [[round(a, 5), round(b, 5)] for a, b in z["poly"]]}
                          for z in geo.SHIPPING_ZONES],
        # lane-comparison corridors per course (sweep.LANES)
        "lanes": {str(c): {leg: {name: [[round(a, 5), round(b, 5)] for a, b in wps]
                                 for name, wps in legs.items() if name != "rhumb line"}
                           for leg, legs in v.items()}
                  for c, v in sweep.LANES.items()},
        "laneHalfWidthM": sweep.LANE_HALF_WIDTH_M,
        "shoal": [[round(a, 5), round(b, 5)] for a, b in geo.BERKELEY_SHOAL],
        "stations": [{"id": s[0], "name": s[1], "lat": s[2], "lon": s[3],
                      "dirFlood": s[4], "dirEbb": s[5]}
                     for s in currents.STATIONS],
        "derivedZones": [{"name": z[0], "lon": z[1], "lat": z[2], "base": z[3],
                          "ampFlood": z[4], "dirFlood": z[5], "ampEbb": z[6],
                          "dirEbb": z[7], "lagMin": z[8]}
                         for z in currents.DERIVED_ZONES],
        "metStations": [{"id": m[0], "name": m[1], "lat": m[2], "lon": m[3]}
                        for m in currents.MET_STATIONS],
        "currentParams": {"idwD0": currents.IDW_D0,
                          "idwPower": currents.IDW_POWER,
                          "shoreReliefL": currents.SHORE_RELIEF_L,
                          "shoreReliefFloor": currents.SHORE_RELIEF_FLOOR},
        "windPoints": [{"name": p[0], "lon": p[1], "lat": p[2],
                        "mult": p[3], "dir": p[4]} for p in wind.WIND_POINTS],
        "windParams": {"idwD0": wind.WIND_IDW_D0, "idwPower": wind.WIND_IDW_POWER,
                       "leePerHeight": wind.LEE_PER_HEIGHT,
                       "leeMax": wind.LEE_MAX_M, "leeDeficit": wind.LEE_DEFICIT,
                       "leeHRef": wind.LEE_H_REF, "leeExp": wind.LEE_EXP},
        "polar": {"tws": polars.TWS_GRID.tolist(), "twa": polars.TWA_GRID.tolist(),
                  "speed": polars.SPEED.tolist(), "name": "Olson 25 (constructed)"},
        "chopK": 0.035,
    }

    # offline fallback: whatever harmonic series we have cached for race day
    date = dt.date.fromisoformat(data["meta"]["date"])
    fb = {}
    for sid, name, lat, lon, dfl, deb in currents.STATIONS:
        try:
            t, v = currents.fetch_noaa_series(date, station=sid)
            fb[sid] = {"t": [round(float(x), 4) for x in t],
                       "v": [round(float(x), 3) for x in v],
                       "dirFlood": dfl, "dirEbb": deb}
        except Exception:
            pass
    data["fallbackSeries"] = fb
    data["fallbackNote"] = ("NOAA harmonic predictions for the race date, "
                            "lst_ldt, + = flood. Used when the live fetch fails.")

    out = HERE / "web" / "model_data.json"
    out.write_text(json.dumps(data, separators=(",", ":")))
    print(f"  web/model_data.json  {out.stat().st_size/1024:.1f} KB  "
          f"({len(data['stations'])} stations, {len(fb)} with cached series, "
          f"{len(data['derivedZones'])} derived zones)")


if __name__ == "__main__":
    main()
