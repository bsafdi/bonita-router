#!/usr/bin/env python3
"""
Fetch NOAA SFBOFS depth-averaged currents and write the JSON the web app reads.

SFBOFS is OPeNDAP/netCDF, which a browser cannot fetch, so this runs on your
laptop race morning.  Forecasts only reach ~48 h ahead.

    pip install xarray netCDF4 numpy
    python prefetch_sfbofs.py --date 2026-09-12 --out sfbofs.json
"""
import argparse, json, datetime as dt
import numpy as np

TEMPLATES = [
    "https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/SFBOFS/MODELS/"
    "{Y}/{M}/{D}/nos.sfbofs.regulargrid.f{f:03d}.{Y}{M}{D}.t{cc}z.nc",
    "https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/SFBOFS/MODELS/"
    "{Y}/{M}/{D}/nos.sfbofs.regulargrid.n{f:03d}.{Y}{M}{D}.t{cc}z.nc",
]
BBOX = dict(lon=(-122.58, -122.29), lat=(37.77, 37.92))
KT = 0.514444


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default=dt.date.today().isoformat())
    ap.add_argument("--cycle", default="03", choices=["03", "09", "15", "21"])
    ap.add_argument("--hours", default="8-21", help="local hour range to keep")
    ap.add_argument("--out", default="sfbofs.json")
    a = ap.parse_args()

    import xarray as xr
    d = dt.date.fromisoformat(a.date)
    Y, M, D = d.strftime("%Y"), d.strftime("%m"), d.strftime("%d")
    h0, h1 = (int(x) for x in a.hours.split("-"))
    frames = []
    for f in range(0, 49):
        ds = None
        for t in TEMPLATES:
            try:
                ds = xr.open_dataset(t.format(Y=Y, M=M, D=D, f=f, cc=a.cycle))
                break
            except Exception:
                pass
        if ds is None:
            continue
        u = np.squeeze(ds[("u_sur" if "u_sur" in ds else "u")].values) / KT
        v = np.squeeze(ds[("v_sur" if "v_sur" in ds else "v")].values) / KT
        lon, lat = ds["lon"].values, ds["lat"].values
        if lon.ndim == 1:
            lon, lat = np.meshgrid(lon, lat)
        m = ((lon >= BBOX["lon"][0]) & (lon <= BBOX["lon"][1]) &
             (lat >= BBOX["lat"][0]) & (lat <= BBOX["lat"][1]) &
             np.isfinite(u) & np.isfinite(v))
        try:
            tval = np.datetime64(ds["time"].values.ravel()[0]).astype("datetime64[m]")
            local = (tval.astype(object) - dt.timedelta(hours=7))
            hour = local.hour + local.minute / 60
        except Exception:
            hour = None
        ds.close()
        if hour is None or not (h0 <= hour <= h1):
            continue
        frames.append({"t": round(float(hour), 3),
                       "lon": [round(float(x), 4) for x in lon[m]],
                       "lat": [round(float(x), 4) for x in lat[m]],
                       "u": [round(float(x), 3) for x in u[m]],
                       "v": [round(float(x), 3) for x in v[m]]})
        print(f"  frame f{f:03d}  local {hour:.2f} h  {int(m.sum())} cells")

    if not frames:
        raise SystemExit("no SFBOFS frames — forecasts only exist ~48 h ahead; "
                         "run this the day before or the morning of the race")
    json.dump({"source": "NOAA SFBOFS", "date": a.date, "cycle": a.cycle,
               "units": "knots", "frames": frames}, open(a.out, "w"),
              separators=(",", ":"))
    print(f"wrote {a.out}: {len(frames)} frames")


if __name__ == "__main__":
    main()
