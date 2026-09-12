"""
Olson 25 velocity prediction.

IMPORTANT: there is no published measured polar for the Olson 25 that I would
trust, so the table below is CONSTRUCTED -- a light-displacement 25-footer that
is hull-speed limited upwind (LWL 21.5 ft -> ~6.2 kt) and planes downwind in a
breeze.  It is shaped to land near the boat's SF Bay PHRF rating (implied ~144
against a real 141-150), and `implied_phrf()` lets you check that.

Shape notes, because these drive the routing answer:
  * The TWA=30 column is NOT zero.  Leaving it at zero made the interpolant
    ramp to nothing between 40 and 30 degrees, which starved every grid course
    inside the tacking cone and made a >2.5 kt foul current look like a wall
    rather than like slow water.  It now sits at ~78% of the 40-degree value,
    which is roughly what the boat does pinching.
  * The planing transition is at 15-16 kt TWS, not 12-13.  A 2500 lb boat with
    a symmetric kite semi-planes at 16 and is fully up by 18-20.
  * Downwind VMG optimum deepens with breeze -- about 150 degrees in the light
    and 160 in 20+ -- which is how a symmetric-spinnaker boat actually runs.
    A fixed 150 at every wind speed was wrong.

If you have real numbers -- a season of tracker data, or a sailmaker's polar --
replace TWS_GRID / TWA_GRID / SPEED and everything downstream just works.
`Polar(scale=...)` is a global speed multiplier for quick calibration.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator

TWS_GRID = np.array([0, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 25, 30], float)
TWA_GRID = np.array([0, 30, 40, 45, 52, 60, 70, 80, 90, 100, 110, 120,
                     130, 140, 150, 160, 170, 180], float)

# rows = TWS_GRID, cols = TWA_GRID, boat speed in knots
SPEED = np.array([
    # 0   30   40   45   52   60   70   80   90  100  110  120  130  140  150  160  170  180
    [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0, 0.0, 2.4, 2.7, 3.0, 3.3, 3.5, 3.6, 3.7, 3.7, 3.6, 3.4, 3.2, 3.0, 2.7, 2.4, 2.1, 1.9],
    [0, 3.1, 4.0, 4.3, 4.5, 4.7, 5.0, 5.1, 5.2, 5.2, 5.1, 4.9, 4.6, 4.3, 3.9, 3.5, 3.1, 2.9],
    [0, 3.7, 4.7, 4.9, 5.1, 5.5, 5.8, 6.0, 6.1, 6.1, 6.0, 5.8, 5.6, 5.3, 4.9, 4.4, 4.0, 3.7],
    [0, 4.1, 5.2, 5.4, 5.6, 6.0, 6.4, 6.7, 6.9, 7.0, 7.0, 6.9, 6.7, 6.4, 5.9, 5.3, 4.8, 4.5],
    [0, 4.3, 5.5, 5.7, 5.9, 6.3, 6.8, 7.2, 7.4, 7.5, 7.4, 7.0, 6.9, 6.7, 6.4, 5.9, 5.4, 5.0],
    [0, 4.4, 5.6, 5.7, 6.0, 6.5, 7.0, 7.4, 7.7, 7.9, 8.0, 8.0, 7.9, 7.6, 7.2, 6.6, 6.0, 5.6],
    [0, 4.5, 5.7, 5.8, 6.1, 6.7, 7.3, 7.8, 8.2, 8.5, 8.7, 8.8, 8.7, 8.4, 8.0, 7.4, 6.7, 6.2],
    [0, 4.5, 5.8, 5.9, 6.2, 6.9, 7.5, 8.1, 8.6, 9.1, 9.6, 10.0, 10.0, 9.8, 9.4, 8.8, 8.1, 7.4],
    [0, 4.5, 5.8, 5.9, 6.2, 7.0, 7.6, 8.3, 8.9, 9.5, 10.2, 10.7, 10.8, 10.6, 10.2, 9.6, 8.9, 8.2],
    [0, 4.4, 5.7, 5.9, 6.2, 7.0, 7.7, 8.4, 9.1, 9.8, 10.6, 11.2, 11.3, 11.1, 10.7, 10.1, 9.4, 8.7],
    [0, 4.3, 5.6, 5.8, 6.1, 7.0, 7.7, 8.5, 9.3, 10.1, 11.0, 11.7, 11.9, 11.7, 11.3, 10.7, 10.0, 9.2],
    [0, 4.1, 5.4, 5.6, 5.9, 6.8, 7.6, 8.4, 9.2, 10.0, 10.9, 11.6, 11.8, 11.6, 11.2, 10.6, 9.9, 9.1],
])


class Polar:
    def __init__(self, tws=TWS_GRID, twa=TWA_GRID, speed=SPEED, scale=1.0,
                 name="Olson 25 (constructed)"):
        self.name = name
        self.tws, self.twa = np.asarray(tws), np.asarray(twa)
        self.speed = np.asarray(speed) * float(scale)
        self.scale = float(scale)
        self._f = RegularGridInterpolator(
            (self.tws, self.twa), self.speed,
            bounds_error=False, fill_value=None)

    def __call__(self, twa_deg, tws_kt):
        """Boat speed through the water.  twa_deg may be any shape; sign ignored."""
        a = np.abs(((np.asarray(twa_deg, float) + 180.0) % 360.0) - 180.0)
        w = np.clip(np.asarray(tws_kt, float), self.tws[0], self.tws[-1])
        shape = np.broadcast_shapes(a.shape, w.shape)
        a = np.broadcast_to(a, shape).ravel()
        w = np.broadcast_to(w, shape).ravel()
        out = np.maximum(self._f(np.stack([w, a], axis=-1)), 0.0)
        return out.reshape(shape) if shape else float(out[0])

    # -- derived quantities -------------------------------------------------

    def best_vmg(self, tws, upwind=True, n=181):
        a = np.linspace(0, 180, n)
        v = self(a, tws)
        vmg = v * np.cos(np.deg2rad(a))
        if upwind:
            i = int(np.argmax(vmg))
        else:
            i = int(np.argmin(vmg))
        return float(a[i]), float(v[i]), float(abs(vmg[i]))

    # Seconds per mile a notional PHRF-0 boat takes round the same triangle.
    # Crude, but it turns the polar into a number you can compare to a rating
    # certificate.  Adjust if you want the check to bite harder.
    SCRATCH_SEC_PER_MILE = 520.0

    def implied_phrf(self, ref_tws=12.0):
        """
        Rough PHRF-equivalent from time round an equilateral triangle (one mile
        each of beat, beam reach and run) in `ref_tws`.  A sanity check that the
        table is not wildly off -- not a rating.
        """
        _, _, vmg_up = self.best_vmg(ref_tws, upwind=True)
        _, _, vmg_dn = self.best_vmg(ref_tws, upwind=False)
        reach = float(self(90.0, ref_tws))
        hours = 1 / vmg_up + 1 / reach + 1 / vmg_dn     # 3 miles sailed made-good
        sec_per_mile = hours * 3600.0 / 3.0
        return sec_per_mile - self.SCRATCH_SEC_PER_MILE


# --------------------------------------------------------------------------
# Wind-against-tide chop
# --------------------------------------------------------------------------

def chop_factor(twd_from, tws, cur_u, cur_v, cog_deg, k=0.045,
                twa_cut=95.0, beat_ref=40.0, run_floor=0.25,
                wind_chop_k=0.05, wind_chop_tws=18.0):
    """
    Speed multiplier for the sea state a light 25-footer actually meets.

    Two terms:

      wind-against-tide  the short steep stuff.  `k` is the fractional loss per
                         knot of fully opposed current AT THE BEAT ANGLE -- the
                         weighting is 1.0 for TWA <= `beat_ref` and falls to
                         zero at `twa_cut`, so 3 kt dead against in 20 kt costs
                         ~18% beating, which is what the Gate ebb does to this
                         boat.  A run through the same chop still loses
                         `run_floor` of that, because the bow stops and the kite
                         collapses; charging it nothing was wrong.
      wind-only chop     the Bay makes its own slop above ~18 kt regardless of
                         the tide.  Small, but it is not zero on a flood day.

    twd_from : radians, direction the wind comes FROM
    cur_u/v  : current in knots, east/north (direction it sets toward)
    cog_deg  : course over ground, degrees true
    """
    wx, wy = -np.sin(twd_from), -np.cos(twd_from)      # unit vector wind blows TOWARD
    cs = np.hypot(cur_u, cur_v)
    with np.errstate(invalid="ignore", divide="ignore"):
        cx = np.where(cs > 1e-6, cur_u / np.maximum(cs, 1e-9), 0.0)
        cy = np.where(cs > 1e-6, cur_v / np.maximum(cs, 1e-9), 0.0)
    opposition = np.clip(-(cx * wx + cy * wy), 0.0, 1.0)   # 1 = dead against
    twa = np.abs(((cog_deg - np.rad2deg(twd_from) + 180.0) % 360.0) - 180.0)

    # 1.0 at and below the beat angle, tapering to zero at twa_cut, never
    # below run_floor -- chop costs something on every point of sail.
    w = np.clip((twa_cut - twa) / max(twa_cut - beat_ref, 1e-6), 0.0, 1.0)
    w = np.maximum(w, run_floor)

    gain = np.clip(tws / 15.0, 0.4, 1.6)               # steeper seas in more breeze
    loss = k * gain * w * opposition * cs
    loss = loss + wind_chop_k * w * np.clip(
        (tws - wind_chop_tws) / 10.0, 0.0, 1.0)
    return np.clip(1.0 - loss, 0.55, 1.0)
