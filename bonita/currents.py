"""
Tidal current field for the Central Bay.

Two back-ends, selected by `source=`:

  "station"  (default, always available)
        A two-tier harmonic model.

          MEASURED  the eight real NOAA current-prediction stations that lie
                    inside the race area.  Each brings its own harmonic series
                    and its own published mean flood/ebb set, so nothing about
                    their amplitude or phase is assumed.
          DERIVED   water NOAA does not instrument -- Point Bonita, the Marin
                    shore, Raccoon Strait, the Berkeley Circle, and the two
                    back-eddies that decide this race.  Each scales a named
                    base station's signed series by a hand amplitude and set.

        The derived tier is the editorial part of the model and the first
        thing to tune.  Both tiers blend by inverse-distance weighting over
        distance THROUGH WATER, so Sausalito never borrows from the Cityfront
        across the headlands.

  "sfbofs"   (high fidelity, only useful within ~48 h of the race)
        NOAA's San Francisco Bay Operational Forecast System: a real 2-D
        hydrodynamic model.  Loader below; see `load_sfbofs`.

This module is the single source of truth for the station and zone tables;
`build_web.py` exports them to the browser app so the two implementations
cannot drift apart.
"""

from __future__ import annotations

import json
import os
import datetime as dt

import numpy as np

from .geo import ll_to_xy, water_distance_field

KT = 0.514444  # m/s per knot

CACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "cache")

# Fallback reference station for the legacy single-station helpers.
REF_STATION = "SFB1201"


# --------------------------------------------------------------------------
# MEASURED tier: NOAA current-prediction stations inside the race area.
# Positions and mean sets are NOAA's own (CO-OPS mdapi / datagetter).
#   id, name, lat, lon, mean flood set, mean ebb set   [degrees true]
# --------------------------------------------------------------------------

STATIONS = [
    ("SFB1201", "San Francisco Bay Entrance",     37.8106, -122.5020,  61, 239),
    ("SFB1202", "Golden Gate Bridge",             37.8292, -122.4620,  52, 238),
    ("SFB1203", "Golden Gate Bridge, 0.46 nm E",  37.8201, -122.4730,  69, 257),
    ("SFB1204", "Alcatraz Island, southwest of",  37.8143, -122.4320,  86, 271),
    ("SFB1205", "Pier 35, north of",              37.8142, -122.4070, 106, 300),
    ("SFB1206", "Pier 23",                        37.8053, -122.3973, 143, 323),
    ("SFB1209", "YBI, W. of Midchannel",          37.8100, -122.3831, 143, 331),
    ("SFB1210", "Treasure Island, 0.78 NM NW of", 37.8373, -122.3872, 126, 302),
]

STATION_IDS = [s[0] for s in STATIONS]


# --------------------------------------------------------------------------
# DERIVED tier  --  EDIT ME
#
# Each zone scales the signed series of `base`.  amp_* is a fraction of that
# station's speed; dir_* is the direction the current SETS TOWARD; lag_min is
# how far this zone lags its base station.
#
# The two rows that matter most are `crissy` and `yellow_bluff`: the Cityfront
# flood back-eddy (setting WEST while the Gate floods east) and the Yellow
# Bluff ebb back-eddy (setting EAST while the Gate ebbs west).  Neither has a
# NOAA station; both decide this race.
# --------------------------------------------------------------------------

DERIVED_ZONES = [
    # name              lon        lat      base       amp_fl dir_fl amp_eb dir_eb lag
    ("bonita_offshore", -122.5500, 37.8050, "SFB1201", 0.75,  80,   0.85, 258,  -5),
    ("bonita_cove",     -122.5250, 37.8210, "SFB1201", 0.55,  55,   0.60, 235,   0),
    ("mile_rocks",      -122.5100, 37.7920, "SFB1201", 0.85,  55,   0.90, 240,   0),
    ("pt_diablo",       -122.5010, 37.8260, "SFB1201", 0.90,  80,   1.00, 250,   2),
    ("gate_north",      -122.4790, 37.8290, "SFB1203", 1.00,  85,   1.05, 265,   0),
    ("gate_south",      -122.4770, 37.8110, "SFB1203", 1.00,  75,   1.10, 262,   0),
    ("yellow_bluff",    -122.4770, 37.8375, "SFB1203", 0.70,  60,   0.40,  90,   8),
    ("sausalito",       -122.4790, 37.8480, "SFB1203", 0.50,  40,   0.30, 110,  12),
    ("crissy",          -122.4650, 37.8045, "SFB1203", 0.40, 250,   0.80, 265,   5),
    ("fort_mason",      -122.4310, 37.8062, "SFB1204", 0.40, 245,   0.85, 270,   6),
    ("harding_rock",    -122.4450, 37.8250, "SFB1203", 0.85,  80,   0.95, 265,   4),
    ("alcatraz_n",      -122.4230, 37.8330, "SFB1204", 0.90,  75,   0.90, 265,   2),
    ("raccoon",         -122.4400, 37.8620, "SFB1204", 0.85,  90,   0.90, 270,  10),
    ("pt_blunt",        -122.4180, 37.8500, "SFB1204", 0.75,  60,   0.80, 250,   8),
    ("angel_ne",        -122.4200, 37.8720, "SFB1210", 0.60,  40,   0.65, 225,   5),
    ("southampton",     -122.4100, 37.8800, "SFB1210", 0.70,  40,   0.75, 220,  10),
    ("berkeley_circle", -122.3720, 37.8600, "SFB1210", 0.60,  60,   0.65, 240,  15),
    ("berkeley_flats",  -122.3350, 37.8650, "SFB1210", 0.35,  70,   0.40, 250,  25),
    ("richmond_s",      -122.3900, 37.8950, "SFB1210", 0.55,  30,   0.60, 210,  15),
]


# Live observed-wind anemometers (CO-OPS met stations, 6-minute).  The browser
# app uses these to anchor the wind field; exported by build_web.py.
MET_STATIONS = [
    ("9414290", "San Francisco (Fort Point)", 37.8063, -122.4659),
    ("9414863", "Richmond",                   37.9283, -122.4000),
    ("9414750", "Alameda",                    37.7720, -122.3003),
]


# Interpolation length scale for blending zones (metres).  Smaller = sharper
# eddies and more local structure; larger = smoother field.
IDW_D0 = 350.0
IDW_POWER = 2.0

# Extra shoreline relief applied on top of the zone table: current magnitude is
# multiplied by (floor + (1-floor)*tanh(d/L)) where d is distance to land.
SHORE_RELIEF_L = 220.0
SHORE_RELIEF_FLOOR = 0.55


# --------------------------------------------------------------------------
# Reference time series
# --------------------------------------------------------------------------

def fetch_noaa_series(date, station=REF_STATION, interval=30, use_cache=True):
    """
    Hourly/half-hourly signed current predictions (knots) for one day.

    Returns (t_hours, v_knots) with t in local hours from midnight.
    Cached to data/cache so the model runs offline once fetched.
    """
    datestr = date.strftime("%Y%m%d")
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, f"currents_{station}_{datestr}_{interval}.json")
    if use_cache:
        import glob
        hits = ([path] if os.path.exists(path) else
                sorted(glob.glob(os.path.join(
                    CACHE_DIR, f"currents_{station}_{datestr}_*.json"))))
        if hits:
            with open(hits[0]) as f:
                d = json.load(f)
            got = os.path.basename(hits[0]).rsplit("_", 1)[1].split(".")[0]
            if got != str(interval):
                print(f"  [currents] {station}: cache holds {got}-minute data, "
                      f"{interval} requested (delete the file to refetch)")
            return np.array(d["t"]), np.array(d["v"])

    import requests
    url = ("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"
           "?product=currents_predictions&application=bonita_router"
           f"&begin_date={datestr}&end_date={datestr}&station={station}"
           f"&time_zone=lst_ldt&interval={interval}&units=english&format=json")
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    rows = r.json()["current_predictions"]["cp"]
    t, v = [], []
    for row in rows:
        ts = dt.datetime.strptime(row["Time"], "%Y-%m-%d %H:%M")
        t.append(ts.hour + ts.minute / 60.0)
        v.append(float(row["Velocity_Major"]))
    t, v = np.array(t), np.array(v)
    with open(path, "w") as f:
        json.dump({"t": t.tolist(), "v": v.tolist(), "station": station,
                   "date": datestr}, f)
    return t, v


class ReferenceSeries:
    """Signed reference current, interpolated in time and wrapped at day edges."""

    def __init__(self, t_hours, v_knots, scale=1.0):
        order = np.argsort(t_hours)
        self.t = np.asarray(t_hours, float)[order]
        self.v = np.asarray(v_knots, float)[order] * scale

    def __call__(self, t_hours):
        return np.interp(np.mod(t_hours, 24.0), self.t, self.v,
                         period=24.0)

    @classmethod
    def for_date(cls, date, **kw):
        t, v = fetch_noaa_series(date, **kw)
        return cls(t, v)

    @classmethod
    def synthetic(cls, max_flood=3.2, max_ebb=-3.8,
                  t_max_flood=11.45, t_max_ebb=17.6):
        """Fallback two-constituent fit if you have no network and no cache."""
        t = np.linspace(0, 24, 289)
        T = 12.42
        a = 0.5 * (max_flood - max_ebb)
        b = 0.5 * (max_flood + max_ebb)
        v = b + a * np.cos(2 * np.pi * (t - t_max_flood) / T)
        # crude diurnal inequality so one ebb is stronger, as on most Bay days
        v += 0.25 * a * np.cos(2 * np.pi * (t - t_max_ebb) / 24.84)
        return cls(t, v)


# --------------------------------------------------------------------------
# Station-based spatial field
# --------------------------------------------------------------------------

def fetch_all_series(date, station_ids=None, interval=30, use_cache=True,
                    verbose=False):
    """
    Signed current predictions for every measured station.

    Returns {id: ReferenceSeries}.  A station that cannot be fetched and has no
    cache is simply left out -- the field degrades to the stations it has, and
    derived zones whose base is missing drop with it.
    """
    out = {}
    for sid in (station_ids or STATION_IDS):
        try:
            t, v = fetch_noaa_series(date, station=sid, interval=interval,
                                     use_cache=use_cache)
            out[sid] = ReferenceSeries(t, v)
        except Exception as e:                          # noqa: BLE001
            if verbose:
                print(f"  [currents] {sid}: {type(e).__name__}")
    return out


class StationCurrentField:
    """
    Two-tier harmonic current field.

    MEASURED zones are the NOAA stations themselves: amplitude 1.0 of their own
    series, NOAA's own mean flood/ebb sets, zero lag.  DERIVED zones scale a
    named base station.  Everything blends by inverse-distance weighting over
    distance through water.
    """

    def __init__(self, grid, series, stations=STATIONS, derived=DERIVED_ZONES,
                 idw_d0=IDW_D0, idw_power=IDW_POWER,
                 shore_relief_L=SHORE_RELIEF_L,
                 shore_relief_floor=SHORE_RELIEF_FLOOR,
                 strength=1.0, verbose=False):
        self.grid = grid
        self.strength = float(strength)

        # Accept either {id: series} or a single series (legacy one-station use).
        if not isinstance(series, dict):
            series = {REF_STATION: series}
        self.series = series

        zones = []
        for (sid, name, lat, lon, dfl, deb) in stations:
            if sid in series:
                zones.append(dict(name=name, kind="measured", base=sid, lon=lon,
                                  lat=lat, amp_fl=1.0, dir_fl=dfl,
                                  amp_eb=1.0, dir_eb=deb, lag=0))
        for (name, lon, lat, base, afl, dfl, aeb, deb, lag) in derived:
            if base in series:
                zones.append(dict(name=name, kind="derived", base=base, lon=lon,
                                  lat=lat, amp_fl=afl, dir_fl=dfl,
                                  amp_eb=aeb, dir_eb=deb, lag=lag))
        if not zones:
            raise RuntimeError("no current stations available")
        self.zones = zones

        W = []
        for z in zones:
            if verbose:
                print("  water-distance field for", z["name"])
            D = water_distance_field(grid, z["lon"], z["lat"])
            W.append(1.0 / (D + idw_d0) ** idw_power)
        W = np.stack(W)
        W[~np.isfinite(W)] = 0.0
        self.w = W / np.maximum(W.sum(axis=0, keepdims=True), 1e-30)

        self.base = [z["base"] for z in zones]
        self.amp_fl = np.array([z["amp_fl"] for z in zones])
        self.dir_fl = np.deg2rad(np.array([z["dir_fl"] for z in zones], float))
        self.amp_eb = np.array([z["amp_eb"] for z in zones])
        self.dir_eb = np.deg2rad(np.array([z["dir_eb"] for z in zones], float))
        self.lag_h = np.array([z["lag"] for z in zones]) / 60.0

        d = grid.dist_to_land
        self.relief = shore_relief_floor + (1 - shore_relief_floor) * np.tanh(
            d / shore_relief_L)

    @property
    def n_measured(self):
        return sum(1 for z in self.zones if z["kind"] == "measured")

    def at(self, t_hours):
        """Depth-averaged current at time t.  Returns (u, v) in KNOTS, east/north."""
        vz = np.array([self.series[b](t_hours - lag)
                       for b, lag in zip(self.base, self.lag_h)])
        flood = vz >= 0
        mag = np.where(flood, self.amp_fl * vz, self.amp_eb * (-vz))
        ang = np.where(flood, self.dir_fl, self.dir_eb)
        u = np.tensordot(mag * np.sin(ang), self.w, axes=(0, 0))
        v = np.tensordot(mag * np.cos(ang), self.w, axes=(0, 0))
        return (u * self.relief * self.strength,
                v * self.relief * self.strength)

    def sample(self, t_hours, lon, lat):
        u, v = self.at(t_hours)
        r, c = self.grid.nearest_water_cell(lon, lat)
        return dict(kt=float(np.hypot(u[r, c], v[r, c])),
                    set=float(np.degrees(np.arctan2(u[r, c], v[r, c])) % 360))


# --------------------------------------------------------------------------
# SFBOFS  --  high-fidelity gridded option
# --------------------------------------------------------------------------

SFBOFS_TEMPLATES = [
    # regular-grid product (easiest); cycle CC in 03,09,15,21
    "https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/SFBOFS/MODELS/"
    "{Y}/{M}/{D}/nos.sfbofs.regulargrid.f{f:03d}.{Y}{M}{D}.t{cc}z.nc",
    "https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/SFBOFS/MODELS/"
    "{Y}/{M}/{D}/nos.sfbofs.regulargrid.n{f:03d}.{Y}{M}{D}.t{cc}z.nc",
]


def load_sfbofs(date, hours, cycle="03", verbose=True):
    """
    Fetch SFBOFS depth-averaged surface currents for a list of forecast hours.

    Returns a list of (t_hours, lon2d, lat2d, u_kt, v_kt).  Requires xarray +
    netCDF4/pydap and network access to opendap.co-ops.nos.noaa.gov.  Forecasts
    only exist for roughly the next 48 h, so this is a race-week tool.
    """
    import xarray as xr
    Y, M, D = date.strftime("%Y"), date.strftime("%m"), date.strftime("%d")
    out = []
    for h in hours:
        ds = None
        for tmpl in SFBOFS_TEMPLATES:
            url = tmpl.format(Y=Y, M=M, D=D, f=int(h), cc=cycle)
            try:
                ds = xr.open_dataset(url)
                break
            except Exception as e:                     # noqa: BLE001
                if verbose:
                    print("  miss", url.split("/")[-1], type(e).__name__)
        if ds is None:
            continue
        uname = "u_sur" if "u_sur" in ds else ("u" if "u" in ds else None)
        vname = "v_sur" if "v_sur" in ds else ("v" if "v" in ds else None)
        if uname is None:
            ds.close()
            continue
        u = np.squeeze(ds[uname].values) / KT
        v = np.squeeze(ds[vname].values) / KT
        lon = ds["lon"].values
        lat = ds["lat"].values
        if lon.ndim == 1:
            lon, lat = np.meshgrid(lon, lat)
        out.append((float(h), lon, lat, u, v))
        ds.close()
    return out


class GriddedCurrentField:
    """Wraps SFBOFS output onto the routing grid, with linear time interpolation."""

    def __init__(self, grid, frames, fallback=None):
        from scipy.interpolate import griddata
        self.grid = grid
        self.times = []
        self.U = []
        self.V = []
        pts_out = np.column_stack([grid.LON.ravel(), grid.LAT.ravel()])
        for (t, lon, lat, u, v) in frames:
            ok = np.isfinite(u) & np.isfinite(v)
            pts = np.column_stack([lon[ok].ravel(), lat[ok].ravel()])
            ui = griddata(pts, u[ok].ravel(), pts_out, method="linear")
            vi = griddata(pts, v[ok].ravel(), pts_out, method="linear")
            self.times.append(t)
            self.U.append(np.nan_to_num(ui).reshape(grid.LON.shape))
            self.V.append(np.nan_to_num(vi).reshape(grid.LON.shape))
        self.times = np.array(self.times)
        self.fallback = fallback

    def at(self, t_hours):
        if len(self.times) == 0:
            if self.fallback is None:
                raise RuntimeError("no SFBOFS frames and no fallback field")
            return self.fallback.at(t_hours)
        i = np.clip(np.searchsorted(self.times, t_hours) - 1, 0,
                    len(self.times) - 2)
        t0, t1 = self.times[i], self.times[i + 1]
        f = 0.0 if t1 == t0 else np.clip((t_hours - t0) / (t1 - t0), 0, 1)
        return (self.U[i] * (1 - f) + self.U[i + 1] * f,
                self.V[i] * (1 - f) + self.V[i + 1] * f)


def build_current_field(grid, date, source="station", strength=1.0,
                        sfbofs_hours=range(0, 30), sfbofs_cycle="03",
                        verbose=False, **kw):
    series = fetch_all_series(date, verbose=verbose)
    if not series:
        print("  [currents] no NOAA data or cache; falling back to a synthetic "
              "harmonic series on the reference station only")
        series = {REF_STATION: ReferenceSeries.synthetic()}
    station = StationCurrentField(grid, series, strength=strength,
                                  verbose=verbose, **kw)
    if verbose or len(series) < len(STATION_IDS):
        print(f"  [currents] {station.n_measured}/{len(STATION_IDS)} NOAA stations, "
              f"{len(station.zones) - station.n_measured} derived zones")
    if source == "station":
        return station
    frames = load_sfbofs(date, sfbofs_hours, cycle=sfbofs_cycle, verbose=verbose)
    if not frames:
        print("  [currents] SFBOFS unavailable; falling back to station model")
        return station
    return GriddedCurrentField(grid, frames, fallback=station)
