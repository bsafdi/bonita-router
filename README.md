# Big Windward–Leeward Regatta — optimal route calculator

A time-optimal routing model for the BYC Big Windward–Leeward Regatta: start
on the Berkeley Circle, round either the Point Bonita buoy "PB" (course 1) or
YRA 15 "EASOM" off Yellow Bluff (course 2) — the course is signalled at the
11:00 warning — and finish back at the start line before the 19:00 limit.
Built for an Olson 25, but the boat is one editable table.

The point of the thing is that this course is decided by current, not by
boatspeed. You beat 12 nm west into the Gate and run 12 nm back, and the tide
turns underneath you halfway through. So the model puts most of its effort
into a defensible current field and then solves for the fastest path through
it.

```
pip install -r requirements.txt
python run_race.py
```

Runs in about 20 s for a single route, ~5 min with the full sensitivity study.

---

## What it computes

**Geometry** (`bonita/geo.py`) — a hand-digitized Central Bay coastline at
~100–300 m fidelity, on a projected grid (250 m by default). Landmasses carry
an effective height, which is what sets how far their wind shadows reach. The
Berkeley flats are an optional shoal mask.

**SI 10 obstructions** (`geo.OBSTRUCTIONS`) — the six lines the Sailing
Instructions make hard no-go (crossing one means retiring): Berkeley Pier
ruins daymark → shore, Point Blunt buoy "3" → Angel Island, Alcatraz bell buoy
AZ → Alcatraz, the "H" beam → shore, Anita Rock (light + AR buoy) → shore, and
the Golden Gate south tower → Fort Point. Five of the six object positions
are taken from the NOAA ENC (chart 18649/18650/18653 data); the H beam is
uncharted and placed from the StFYC SI text (~200 yd W of the club). Each line
is finished at the closest point of the land polygon and sealed 200 m inland.
Enforcement is exact, not rasterised: every grid edge whose segment crosses a
line is removed from the graph (all 16 move directions), and cells within 70 m
of a line are closed. `tests/test_solver.py` proves every line uncrossable.

**Shipping lanes** (`geo.SHIPPING_ZONES`, SI 9.1 / Inland Rule 9) — the Main
Ship Channel lanes, the two SF Bay lanes either side of Alcatraz, and the
Golden Gate and Central Bay precautionary areas, from the ENC approach cells.
These are *soft*: time inside a lane is charged ×1.10 (Gate precautionary
area ×1.05, Central Bay flag-only) in the search objective, real time is
reported, and every crossing is listed with its clock times so the tactician
knows where a ship can force them off. `--no-lane-cost` keeps the flags only.

**Currents** (`bonita/currents.py`) — two back-ends:

- `station` (default): one NOAA harmonic reference station (`SFB1201`, the Bay
  entrance) supplies the time dependence; an editable table of ~24 zones
  supplies the spatial structure. Each zone carries a *separate* flood and ebb
  amplitude and set, plus a phase lag. That separation is the important part —
  it is what reproduces the **Cityfront flood back-eddy** (setting *west* at
  ~0.5 kt while the Gate floods at 3 kt) and the **Yellow Bluff ebb
  back-eddy** (setting *east* while the Gate ebbs). Zones are blended by
  inverse-distance weighting over *water* distance, so Sausalito never borrows
  from the Cityfront across the headlands.
- `sfbofs`: NOAA's San Francisco Bay Operational Forecast System — a real 2-D
  hydrodynamic model. Forecasts only run ~48 h ahead, so this is the race-week
  option. `--currents sfbofs` falls back to the station model if it can't
  fetch.

**Wind** (`bonita/wind.py`) — two back-ends:

- `param` (default): a diurnal sea-breeze curve (`base → peak → decay`)
  times a static spatial structure. The structure is a table of control points
  (Gate acceleration, the slot, the Cityfront, the Berkeley shore) *plus* a
  directional terrain shadow computed by ray-marching upwind from every cell.
  The shadow is derived, not hand-placed, so the Angel Island, Alcatraz and
  Sausalito lees re-derive themselves if you change the gradient direction.
- `openmeteo`: real gridded forecast, no API key. ~3 km resolution smears the
  Gate acceleration, so the terrain structure is re-applied on top as a
  normalized overlay (`terrain_overlay=False` to disable).

**Boat** (`bonita/polars.py`) — an Olson 25 polar, plus a wind-against-tide
chop penalty that bites upwind and not downwind.

**Router** (`bonita/router.py`) — earliest-arrival Dijkstra over
`(grid cell, direction of arrival)` states. Carrying the arrival direction is
what lets it charge a real tack/gybe penalty, which is what stops a grid router
from producing a physically meaningless sawtooth: with the penalty in, a beat
comes out as a handful of long tacks. Every edge is priced by solving the
current triangle properly — the grid fixes the course over ground, and the
solver searches boat headings for the one whose water velocity plus the tidal
set points along that COG.

**Sweeps** (`bonita/sweep.py`) — lane comparison (each strategy sailed
*optimally within its own corridor*, so it measures the cost of committing to a
side rather than the cost of steering badly), start time, breeze scenario, and
polar scale.

---

## Updating it as the race gets closer

Roughly in order of how much they matter.

**1. Currents — refresh, then upgrade.** The reference series is cached in
`data/cache/`. Delete the file to re-fetch, or just run on a different
`--date`. Within ~48 h of the race, switch to the real model:

```bash
python run_race.py --currents sfbofs
```

**2. Wind — put the forecast in.** Two ways, cheapest first. Set the sea-breeze
knobs by hand from whatever forecast you trust:

```bash
python run_race.py --peak-kt 22 --onset 11:00 --wind-peak-time 16:00
```

Or pull a gridded forecast (only meaningful inside ~3 days):

```bash
python run_race.py --wind openmeteo
```

**3. The polar.** This is the weakest input — see the caveats below. If you get
real numbers, replace the table in `polars.py`. In the meantime,
`--polar-scale` is a global multiplier and `Polar().implied_phrf()` tells you
where the table sits against the rating.

**4. Marks.** `MARKS` in `geo.py`. PB and EASOM are the SI positions. The
start is an estimate on the DOC–FOC line near DOC (SI 5.1 says "~1.5 NM W of
the breakwater, between DOC and FOC"); put the signal boat's real position in
when you see it — everything downstream follows.

Useful flags:

| flag | what it does |
|---|---|
| `--course 1` / `2` / `both` | SI 12.0 course: Point Bonita (1) or EASOM (2); default both |
| `--start 11:05` | start gun (default: SI 3.1 warning 11:00 + 5 min) |
| `--no-obstructions` | ignore the SI 10 lines (comparison only) |
| `--no-lane-cost` | flag shipping-lane crossings without penalising them |
| `--skip-unconstrained` | do not also solve the unconstrained route for comparison |
| `--res 200` | finer grid (slower; 350 for a quick look) |
| `--current-scale 1.15` | scale the whole tidal field — spring/neap fudge |
| `--avoid-shoal` | treat the Berkeley flats as unnavigable |
| `--only route` / `--only sweeps` | skip half the work |

---

## What it gets right, and what to distrust

Ranked by how much I'd trust it.

**Solid.** The routing algorithm and the current-triangle solve are exact given
the fields — no hand-waving there. The *shape* of the tidal argument is
robust: the reference series is real NOAA harmonic data, and the flood/ebb eddy
structure is standard Bay knowledge, not invented.

**Reasonable but tunable.** The current zone table. The amplitudes and sets are
calibrated by eye against NOAA current charts. Individual zones could be off by
20–30%, which moves lane comparisons by a few minutes but does not flip the
big conclusions.

**Distrust.** Three things, in order:

1. **The polar is constructed.** There is no measured Olson 25 polar I'd trust,
   so I built one shaped to land near the boat's PHRF (`implied_phrf()` ≈ 153;
   the SF Bay rating is around 141–150, so it's a touch conservative). The
   sensitivity is steep — a 4% speed change moves the implied rating by ~30
   s/mile and the elapsed time by ~10 minutes. **Absolute elapsed times are
   only as good as this table.** Relative comparisons — which lane, how much a
   later start helps — are far more robust, because the polar error largely
   cancels.

2. **The coastline is hand-digitized**, not surveyed. Fine at 250 m, wrong if
   you zoom in on shoreline-hugging tactics. If you care about the last 100 m
   along the Sausalito or Cityfront shore, replace `LAND_POLYGONS` with real
   ENC or OSM data.

3. **The parameterized wind field is a caricature.** The Gate acceleration and
   the island lees are real effects at roughly the right magnitude, but the
   numbers are judgment, not measurement. It also has no shifts, no gusts, and
   no puff structure — it is a smooth deterministic field, so the router never
   has to make a decision under uncertainty, which is most of what real
   tactics is.

**Structurally absent.** No fleet, so no starting-line effects, no dirty air,
no covering, no right-of-way. No waves as a separate state (only the chop
penalty). The tidal current is depth-averaged, so no vertical shear. And the
optimizer is deterministic: it knows the future exactly, which is worth
something real on a course where the tide turns mid-race. Treat its route as
the best case a boat with perfect foreknowledge could sail, and its elapsed
time as an optimistic bound.

---

## Layout

```
run_race.py            CLI
bonita/geo.py          coastline, grid, marks, land heights
bonita/currents.py     NOAA station model + zone table + SFBOFS loader
bonita/wind.py         sea-breeze model + terrain shadow + Open-Meteo loader
bonita/polars.py       Olson 25 polar, wind-against-tide chop
bonita/router.py       Dijkstra router, current-triangle solve
bonita/sweep.py        lane comparison and sensitivity studies
bonita/plotting.py     charts
data/cache/            fetched NOAA / forecast data (git-ignorable)
out/                   figures and run logs
```
