"""Charts: route over the current field, route over the wind field, leg tables."""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402

from .geo import MARKS, ll_to_xy, SHIPPING_ZONES, obstruction_lines_ll  # noqa: E402

ASPECT = 1.0 / np.cos(np.deg2rad(37.85))


def hhmm(t):
    t = float(t) % 24.0
    h = int(t)
    m = int(round((t - h) * 60))
    if m == 60:
        h, m = h + 1, 0
    return f"{h:02d}:{m:02d}"


def _basemap(ax, grid):
    ax.pcolormesh(grid.LON, grid.LAT, grid.land, cmap="Greys",
                  vmin=-0.6, vmax=2.2, shading="auto", zorder=0)
    ax.set_aspect(ASPECT)
    ax.set_xlim(-122.565, -122.305)
    ax.set_ylim(37.782, 37.905)
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")


def plot_route(grid, current_field, wind_field, legs, t_snapshots=None,
               title="", path="route.png", show="current"):
    """
    legs: list of RouteResult.  show = "current" | "wind".
    One panel per snapshot time, route drawn on all of them, with the part of
    the route sailed before that time highlighted.
    """
    ts = t_snapshots or [legs[0].t_start + 0.5 * legs[0].elapsed_h]
    n = len(ts)
    fig, axes = plt.subplots(1, n, figsize=(8.0 * n, 6.4), squeeze=False)
    for ax, t in zip(axes[0], ts):
        _basemap(ax, grid)
        if show == "current":
            u, v = current_field.at(t)
            sp = np.where(grid.water, np.hypot(u, v), np.nan)
            m = ax.pcolormesh(grid.LON, grid.LAT, sp, cmap="magma",
                              vmin=0, vmax=3.8, shading="auto", zorder=1)
            lab, s = "current (kt)", 3
            U = np.where(grid.water, u, np.nan)
            V = np.where(grid.water, v, np.nan)
            qs = 42
        else:
            sp_, dr = wind_field.at(t)
            sp = np.where(grid.water, sp_, np.nan)
            m = ax.pcolormesh(grid.LON, grid.LAT, sp, cmap="viridis",
                              vmin=4, vmax=28, shading="auto", zorder=1)
            lab, s = "TWS (kt)", 3
            U = np.where(grid.water, -np.sin(dr) * sp_, np.nan)
            V = np.where(grid.water, -np.cos(dr) * sp_, np.nan)
            qs = 430
        ax.quiver(grid.LON[::s, ::s], grid.LAT[::s, ::s],
                  U[::s, ::s], V[::s, ::s], scale=qs, width=0.0022,
                  color="w", alpha=0.65, zorder=2)
        plt.colorbar(m, ax=ax, label=lab, fraction=0.035, pad=0.02)

        for leg, col in zip(legs, ["#00e5ff", "#ff4d6d"]):
            if leg is None or not leg.ok:
                continue
            ax.plot(leg.lon, leg.lat, "-", color=col, lw=1.7, zorder=3)
            done = np.array(leg.path_t) <= t
            if done.any():
                ax.plot(leg.lon[done], leg.lat[done], "-", color="w", lw=3.2,
                        alpha=0.85, zorder=4)
                ax.plot(leg.lon[done][-1], leg.lat[done][-1], "o", color="w",
                        ms=7, mec="k", zorder=6)
        for nme, (lo, la) in MARKS.items():
            ax.plot(lo, la, "*", color="yellow", ms=17, mec="k", zorder=7)
        # SI 9.1 shipping zones (soft) and SI 10 obstruction lines (hard)
        for z in SHIPPING_ZONES:
            P = np.array(z["poly"] + [z["poly"][0]])
            ax.plot(P[:, 0], P[:, 1], "--", color="#c8c8ff", lw=0.8,
                    alpha=0.7, zorder=2.5)
        for _key, pts in obstruction_lines_ll():
            P = np.array(pts)
            ax.plot(P[:, 0], P[:, 1], "-", color="#ff2020", lw=2.2, zorder=5)
        ax.set_title(f"{hhmm(t)}", fontsize=12)
    fig.suptitle(title, fontsize=13)
    fig.tight_layout()
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return path


def leg_table(router, leg, every_min=20):
    """Position / speed / wind / current every `every_min` along a leg."""
    rows = []
    t_prev = None
    for (r, c), t in zip(leg.path_rc, leg.path_t):
        if t_prev is not None and (t - t_prev) * 60 < every_min:
            continue
        t_prev = t
        _b, (tws, twd, cu, cv) = router.fields(t)
        rows.append(dict(
            time=hhmm(t),
            lon=float(router.g.LON[r, c]), lat=float(router.g.LAT[r, c]),
            tws=float(tws[r, c]),
            twd=float(np.rad2deg(twd[r, c]) % 360),
            cur_kt=float(np.hypot(cu[r, c], cv[r, c])),
            cur_set=float(np.rad2deg(np.arctan2(cu[r, c], cv[r, c])) % 360),
        ))
    return rows


def track_length_nm(leg):
    x, y = ll_to_xy(leg.lon, leg.lat)
    return float(np.sum(np.hypot(np.diff(x), np.diff(y))) / 1852.0)


def plot_sweep(df, path="sweep.png", title="Start-time sensitivity"):
    """df: list of dicts with keys start, total_h, and optionally scenario."""
    import collections
    by = collections.defaultdict(list)
    for row in df:
        by[row.get("scenario", "base")].append(row)
    fig, ax = plt.subplots(figsize=(9, 5.2))
    for name, rows in by.items():
        rows = sorted(rows, key=lambda r: r["start"])
        ax.plot([r["start"] for r in rows],
                [r["total_h"] for r in rows], "o-", label=name)
    ax.set_xlabel("start time (hours, local)")
    ax.set_ylabel("elapsed time (hours)")
    ax.grid(alpha=0.3)
    ax.legend()
    ax.set_title(title)
    xt = ax.get_xticks()
    ax.set_xticks(xt)
    ax.set_xticklabels([hhmm(v) for v in xt])
    ax.set_xlim(min(r["start"] for r in df) - 0.1,
                max(r["start"] for r in df) + 0.1)
    fig.tight_layout()
    fig.savefig(path, dpi=115)
    plt.close(fig)
    return path
