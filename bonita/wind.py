"""
Wind field for the Central Bay.

Two back-ends, selected by `source=`:

  "param"   (default)
        A parameterised sea breeze.  One diurnal strength curve
        (base -> peak -> decay) multiplied by a static spatial structure:
          * a table of control points carrying a speed multiplier and a
            direction offset (the Gate acceleration, the slot, the Cityfront),
          * a *directional* terrain shadow computed by ray-marching upwind
            from every cell, which produces the Angel Island / Alcatraz /
            Treasure Island / Sausalito lees automatically and re-derives them
            if you change the gradient direction.
        Everything is a knob, so this is the back-end for scenario sweeps.

  "openmeteo"  (race week)
        Real gridded forecast (HRRR / ECMWF via api.open-meteo.com), no API
        key.  ~3 km resolution smears the Gate acceleration, so by default the
        terrain structure above is re-applied on top as a normalised overlay.
"""

from __future__ import annotations

import json
import os
import datetime as dt

import numpy as np

from .geo import ll_to_xy

CACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "cache")


# --------------------------------------------------------------------------
# Diurnal strength curve  --  EDIT ME
# --------------------------------------------------------------------------

class SeaBreeze:
    """
    Domain-mean true wind speed vs time of day, in knots.

    base_kt   morning gradient breeze
    peak_kt   afternoon maximum (this is the number you tune from the forecast)
    t_onset   hour the breeze starts filling
    t_peak    hour it reaches peak
    t_hold    hour it starts dying
    t_off     hour it is back to base
    """

    def __init__(self, base_kt=8.0, peak_kt=19.0, t_onset=10.0, t_peak=15.5,
                 t_hold=18.0, t_off=21.0):
        self.base_kt, self.peak_kt = float(base_kt), float(peak_kt)
        self.t_onset, self.t_peak = float(t_onset), float(t_peak)
        self.t_hold, self.t_off = float(t_hold), float(t_off)

    def __call__(self, t):
        t = np.asarray(t, dtype=float)
        s = np.full(t.shape, self.base_kt)
        # raised-cosine build
        m = (t >= self.t_onset) & (t < self.t_peak)
        f = (t - self.t_onset) / max(self.t_peak - self.t_onset, 1e-6)
        s = np.where(m, self.base_kt + (self.peak_kt - self.base_kt) *
                     0.5 * (1 - np.cos(np.pi * np.clip(f, 0, 1))), s)
        s = np.where((t >= self.t_peak) & (t < self.t_hold), self.peak_kt, s)
        m = (t >= self.t_hold) & (t < self.t_off)
        f = (t - self.t_hold) / max(self.t_off - self.t_hold, 1e-6)
        s = np.where(m, self.peak_kt + (self.base_kt - self.peak_kt) *
                     0.5 * (1 - np.cos(np.pi * np.clip(f, 0, 1))), s)
        s = np.where(t >= self.t_off, self.base_kt, s)
        return s if s.shape else float(s)


# --------------------------------------------------------------------------
# Spatial structure  --  EDIT ME
#
# mult    : speed multiplier relative to the domain mean
# dir_deg : direction the wind is coming FROM, degrees true
# --------------------------------------------------------------------------

WIND_POINTS = [
    # name                  lon        lat     mult  dir
    ("ocean_west",      -122.5650, 37.8000,   1.00, 282),
    ("bonita",          -122.5350, 37.8150,   1.02, 278),
    ("headland_lee",    -122.4980, 37.8265,   0.72, 252),   # Kirby Cove / Pt Diablo: light and shifty
    ("gate_throat",     -122.4783, 37.8199,   1.30, 270),   # gap acceleration
    ("gate_south",      -122.4770, 37.8105,   1.26, 268),
    ("gate_north",      -122.4790, 37.8285,   1.18, 272),
    ("crissy",          -122.4650, 37.8050,   1.18, 272),
    ("fort_mason",      -122.4310, 37.8065,   0.92, 256),
    ("aquatic_park",    -122.4225, 37.8062,   0.74, 250),
    ("harding_rock",    -122.4450, 37.8250,   1.20, 265),
    ("alcatraz_w",      -122.4290, 37.8265,   1.14, 262),
    ("yellow_bluff",    -122.4770, 37.8380,   0.88, 258),
    ("sausalito",       -122.4790, 37.8490,   0.58, 240),   # deep lee under the Marin hills
    ("raccoon",         -122.4400, 37.8620,   0.70, 250),
    ("angel_south",     -122.4290, 37.8480,   0.96, 258),
    ("pt_blunt",        -122.4180, 37.8500,   1.06, 255),
    ("southampton",     -122.4100, 37.8800,   1.00, 246),
    ("berkeley_circle", -122.3720, 37.8600,   1.08, 250),   # the slot delivers here
    ("berkeley_inshore",-122.3350, 37.8650,   0.82, 246),   # marina shore effect
    ("richmond_s",      -122.3900, 37.8950,   0.90, 242),
    ("treasure_e",      -122.3560, 37.8220,   0.74, 254),
    ("south_central",   -122.3800, 37.7950,   0.98, 254),
]

WIND_IDW_D0 = 900.0     # metres; larger than for currents -- wind fields are smoother
WIND_IDW_POWER = 2.0

# Directional terrain shadow.  The shadow length scales with the height of the
# upwind land: a 240 m island like Angel throws a lee well over a mile, flat
# Treasure Island throws almost none.
LEE_PER_HEIGHT = 14.0   # shadow length = this many x obstacle height (metres)
LEE_MAX_M = 4000.0      # cap on shadow length
LEE_DEFICIT = 0.55      # max fractional speed loss, for tall land
LEE_H_REF = 150.0       # height at which the deficit saturates
LEE_EXP = 0.60          # recovery exponent with downwind distance


class ParamWindField:
    def __init__(self, grid, breeze=None, points=WIND_POINTS,
                 idw_d0=WIND_IDW_D0, idw_power=WIND_IDW_POWER,
                 lee_per_height=LEE_PER_HEIGHT, lee_max=LEE_MAX_M,
                 lee_deficit=LEE_DEFICIT, lee_h_ref=LEE_H_REF, lee_exp=LEE_EXP,
                 dir_shift_deg=0.0, speed_scale=1.0):
        self.grid = grid
        self.breeze = breeze or SeaBreeze()
        self.speed_scale = float(speed_scale)

        # --- IDW over Euclidean distance (wind crosses land, unlike current)
        W, mult, dirs = [], [], []
        for (_n, lon, lat, m, d) in points:
            px, py = ll_to_xy(lon, lat)
            dist = np.hypot(grid.X - px, grid.Y - py)
            W.append(1.0 / (dist + idw_d0) ** idw_power)
            mult.append(m)
            dirs.append(np.deg2rad(d + dir_shift_deg))
        W = np.stack(W)
        W /= W.sum(axis=0, keepdims=True)
        self.mult = np.tensordot(np.array(mult), W, axes=(0, 0))
        # interpolate direction as a unit vector so it never wraps badly
        cx = np.tensordot(np.cos(np.array(dirs)), W, axes=(0, 0))
        sx = np.tensordot(np.sin(np.array(dirs)), W, axes=(0, 0))
        self.dir_from = np.arctan2(sx, cx)          # radians, FROM

        # --- directional terrain shadow, ray-marched upwind
        self.lee = self._lee_factor(lee_per_height, lee_max, lee_deficit,
                                    lee_h_ref, lee_exp)

        # combined static multiplier, renormalised so the domain mean over
        # water is 1 (so `peak_kt` means what it says)
        combo = self.mult * self.lee
        self.struct = combo / np.mean(combo[grid.water])

    def _lee_factor(self, lee_per_height, lee_max, lee_deficit, lee_h_ref,
                    lee_exp):
        """
        March upwind from every cell.  The first land hit sets the shadow: its
        height fixes both how deep the deficit is and how far downwind it
        reaches, so a tall island shadows far and a flat one barely at all.
        """
        g = self.grid
        ux = np.sin(self.dir_from)      # unit vector toward where wind comes FROM
        uy = np.cos(self.dir_from)
        n_steps = max(int(lee_max / g.res), 1)
        hit_dist = np.full(g.land.shape, np.inf)
        hit_h = np.zeros(g.land.shape)
        x0, y0 = g.x[0], g.y[0]
        for s in range(1, n_steps + 1):
            d = s * g.res
            ci = np.clip(np.round((g.X + ux * d - x0) / g.res).astype(int),
                         0, g.nx - 1)
            ri = np.clip(np.round((g.Y + uy * d - y0) / g.res).astype(int),
                         0, g.ny - 1)
            hit = g.land[ri, ci] & ~np.isfinite(hit_dist)
            hit_dist = np.where(hit, d, hit_dist)
            hit_h = np.where(hit, g.height[ri, ci], hit_h)

        reach = np.minimum(lee_per_height * hit_h, lee_max)
        deficit = lee_deficit * np.clip(hit_h / lee_h_ref, 0, 1)
        with np.errstate(invalid="ignore", divide="ignore"):
            frac = np.clip(hit_dist / np.maximum(reach, 1.0), 0, 1)
        f = np.where(np.isfinite(hit_dist) & (reach > 0),
                     1.0 - deficit * (1.0 - frac ** lee_exp),
                     1.0)
        return f

    def at(self, t_hours):
        """Returns (tws_kt, twd_rad_from) arrays on the grid."""
        s = float(self.breeze(t_hours)) * self.speed_scale
        return s * self.struct, self.dir_from


# --------------------------------------------------------------------------
# Gridded forecast (Open-Meteo)
# --------------------------------------------------------------------------

def fetch_openmeteo(date, lats, lons, model="best_match", use_cache=True):
    """Hourly 10 m wind for a set of points.  Returns dict with arrays."""
    datestr = date.strftime("%Y-%m-%d")
    os.makedirs(CACHE_DIR, exist_ok=True)
    key = f"wind_{model}_{datestr}_{len(lats)}pts.json"
    path = os.path.join(CACHE_DIR, key)
    if use_cache and os.path.exists(path):
        return json.load(open(path))

    import requests
    url = ("https://api.open-meteo.com/v1/forecast"
           f"?latitude={','.join(f'{v:.4f}' for v in lats)}"
           f"&longitude={','.join(f'{v:.4f}' for v in lons)}"
           "&hourly=wind_speed_10m,wind_direction_10m,wind_gusts_10m"
           f"&wind_speed_unit=kn&timezone=America%2FLos_Angeles"
           f"&start_date={datestr}&end_date={datestr}&models={model}")
    r = requests.get(url, timeout=90)
    r.raise_for_status()
    data = r.json()
    if isinstance(data, dict):
        data = [data]
    out = {"lat": list(lats), "lon": list(lons),
           "time": [t[11:16] for t in data[0]["hourly"]["time"]],
           "speed": [d["hourly"]["wind_speed_10m"] for d in data],
           "dir": [d["hourly"]["wind_direction_10m"] for d in data],
           "gust": [d["hourly"].get("wind_gusts_10m") for d in data]}
    json.dump(out, open(path, "w"))
    return out


class GriddedWindField:
    """Open-Meteo points -> smooth field on the routing grid, optional terrain overlay."""

    def __init__(self, grid, data, terrain_overlay=None, speed_scale=1.0,
                 idw_d0=3000.0):
        self.grid = grid
        self.speed_scale = float(speed_scale)
        self.overlay = terrain_overlay
        px, py = ll_to_xy(np.array(data["lon"]), np.array(data["lat"]))
        W = []
        for i in range(len(px)):
            dist = np.hypot(grid.X - px[i], grid.Y - py[i])
            W.append(1.0 / (dist + idw_d0) ** 2)
        W = np.stack(W)
        self.w = W / W.sum(axis=0, keepdims=True)
        self.t = np.array([int(s[:2]) + int(s[3:]) / 60 for s in data["time"]])
        self.spd = np.array(data["speed"])      # (npts, ntime)
        self.dirdeg = np.array(data["dir"])

    def at(self, t_hours):
        i = int(np.clip(np.searchsorted(self.t, t_hours) - 1, 0, len(self.t) - 2))
        f = np.clip((t_hours - self.t[i]) / (self.t[i + 1] - self.t[i]), 0, 1)
        sp = self.spd[:, i] * (1 - f) + self.spd[:, i + 1] * f
        a0 = np.deg2rad(self.dirdeg[:, i])
        a1 = np.deg2rad(self.dirdeg[:, i + 1])
        cx = np.cos(a0) * (1 - f) + np.cos(a1) * f
        sx = np.sin(a0) * (1 - f) + np.sin(a1) * f
        speed = np.tensordot(sp, self.w, axes=(0, 0))
        cxg = np.tensordot(cx, self.w, axes=(0, 0))
        sxg = np.tensordot(sx, self.w, axes=(0, 0))
        dirf = np.arctan2(sxg, cxg)
        if self.overlay is not None:
            # The forecast grid already resolves part of the terrain effect at
            # ~3 km, so applying the full structure on top double counts it.
            # Take the square root: the overlay supplies the sub-grid half.
            speed = speed * np.sqrt(np.maximum(self.overlay.struct, 1e-6))
            # keep the forecast's mean direction but adopt local terrain veer
            dirf = dirf + (self.overlay.dir_from -
                           np.mean(self.overlay.dir_from[self.grid.water]))
        return speed * self.speed_scale, dirf


def build_wind_field(grid, date, source="param", breeze=None,
                     terrain_overlay=True, speed_scale=1.0,
                     n_forecast_pts=12, model="best_match", **kw):
    param = ParamWindField(grid, breeze=breeze, speed_scale=speed_scale, **kw)
    if source == "param":
        return param
    # a coarse scatter of forecast points across the domain water
    lat = np.linspace(37.79, 37.90, 4)
    lon = np.linspace(-122.55, -122.32, 5)
    LA, LO = np.meshgrid(lat, lon)
    pts = [(a, o) for a, o in zip(LA.ravel(), LO.ravel())]
    pts = pts[:n_forecast_pts] if n_forecast_pts else pts
    try:
        data = fetch_openmeteo(date, [p[0] for p in pts], [p[1] for p in pts],
                               model=model)
    except Exception as e:                              # noqa: BLE001
        print(f"  [wind] forecast fetch failed ({type(e).__name__}); "
              f"using parameterised field")
        return param
    return GriddedWindField(grid, data,
                            terrain_overlay=param if terrain_overlay else None,
                            speed_scale=speed_scale)
