# Bonita Router — working notes for Claude Code

Race: **BYC Big Windward–Leeward Regatta, Saturday 12 September 2026.**
First warning **1100**, gun modelled at **11:05**. Time limit **1900** (SI 6.2).
Boat: **Olson 25**, SF Bay PHRF ~141–150.

## Run it

```bash
pip install -r requirements.txt
python run_race.py                 # both courses, 11:05, race-day tide
python run_race.py --course 1      # Point Bonita only
python build_web.py                # rebuild dist/ and web/artifact.html
python -m tests.test_solver        # must stay green
```

## Architecture, and the one rule

**`bonita/` is the single source of truth for every table.** `export_model_data.py`
copies them into `web/model_data.json`; `build_web.py` calls it and inlines the
result. Never hand-edit anything under `web/model_data.json`, and never add a
constant to the JavaScript that does not exist in Python — the two
implementations are verified numerically against each other and must not drift.

| file | holds |
|---|---|
| `bonita/geo.py` | coastline, land heights, grid, `MARKS`, `COURSES`, `OBSTRUCTIONS`, `SHIPPING_ZONES` |
| `bonita/currents.py` | 8 NOAA stations + 19 derived zones, SFBOFS loader |
| `bonita/wind.py` | sea-breeze curve, control points, terrain-shadow ray march |
| `bonita/polars.py` | Olson 25 polar, chop model |
| `bonita/router.py` | Dijkstra, current-triangle edge pricing, obstruction edges |
| `web/model.js` | browser port — mirrors the above line for line |
| `web/app.body.html` | the on-water UI |

## Things that will bite you

- **SI 10 obstructions are hard.** Crossing one means retiring. They are
  enforced by removing every grid *edge* whose segment crosses the line, not by
  masking cells — a 2:1 knight move would otherwise hop a one-cell strip. If
  you touch grid connectivity, re-run the obstruction tests.
- **Water-distance Dijkstra must use Float64.** A Float32 store against
  float64 heap keys silently truncated the search to 20% of the Bay once.
- **`grid.land` is the polygon land only.** The shore buffer narrows `water`,
  not `land` — the wind lee ray-march reads `land` and its height, and folding
  the buffer in erases every wind shadow.
- **CO-OPS `units=english` is already knots.** Do not convert.
- **Edges longer than 3 time bins are refused.** Fields are frozen at the
  departure bin, so a long edge is stale water; the router will otherwise
  "sail" a cell at 0.2 kt for 100 minutes through a flood that has eased.
- Marks are routed *to*, not rounded to port. SI 12.3 says port roundings.

## Known-weak inputs, in order

1. The Olson 25 polar is constructed, not measured (implied PHRF 142). Absolute
   times ±10 min; lane comparisons are far more robust.
2. NOAA's Golden Gate stations may report a depth bin that understates the
   surface ebb by ~1 kt. Would make the run home slower than modelled.
3. The Cityfront flood eddy is one grid cell wide and should be a band of
   4–6 seeds. The Marin-shore ebb zones are editorial, not measured.
4. The "H" beam (SI 10.1.d) position is estimated ±100 m from text; the other
   five obstruction objects are from the NOAA ENC.
5. Start line is an estimate on the YRA DOC–FOC line. Replace with the signal
   boat's real position when you see it: `MARKS["start"]` in `geo.py`.

## Deploying the on-water app

`dist/` is a static site. Push it to GitHub Pages (see `dist/DEPLOY.md`, or use
the included workflow). It fetches NOAA currents, NOAA anemometers and
Open-Meteo directly from the phone — all three allow cross-origin browser
requests, verified. It installs to the home screen and caches for no-signal use.
