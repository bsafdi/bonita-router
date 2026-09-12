/* =====================================================================
   Bonita router -- browser port of the Python model.

   Same physics, same numbers: grid + land mask from the shipped polygons,
   a station-driven current field with flood/ebb eddy zones, a parameterised
   sea breeze with a ray-marched terrain shadow, an Olson 25 polar, and an
   earliest-arrival Dijkstra over (cell, arrival direction) states.

   Everything runs from the ~9 KB model_data.json blob, so the tables stay
   editable rather than being baked into precomputed fields.
   ===================================================================== */
(function (root) {
'use strict';

const D2R = Math.PI / 180, R2D = 180 / Math.PI;

/* ---------------------------------------------------------------- heap */
class MinHeap {
  constructor(cap) {
    this.k = new Float64Array(cap);      // key (time)
    this.v = new Int32Array(cap);        // payload (state index)
    this.n = 0;
    this.cap = cap;
  }
  _grow() {
    const k = new Float64Array(this.cap * 2), v = new Int32Array(this.cap * 2);
    k.set(this.k); v.set(this.v);
    this.k = k; this.v = v; this.cap *= 2;
  }
  push(key, val) {
    if (this.n === this.cap) this._grow();
    let i = this.n++;
    this.k[i] = key; this.v[i] = val;
    while (i > 0) {
      const p = (i - 1) >> 1;
      if (this.k[p] <= this.k[i]) break;
      const tk = this.k[p], tv = this.v[p];
      this.k[p] = this.k[i]; this.v[p] = this.v[i];
      this.k[i] = tk; this.v[i] = tv;
      i = p;
    }
  }
  pop() {
    if (this.n === 0) return null;
    const rk = this.k[0], rv = this.v[0];
    this.n--;
    if (this.n > 0) {
      this.k[0] = this.k[this.n]; this.v[0] = this.v[this.n];
      let i = 0;
      for (;;) {
        const l = 2 * i + 1, r = l + 1;
        let s = i;
        if (l < this.n && this.k[l] < this.k[s]) s = l;
        if (r < this.n && this.k[r] < this.k[s]) s = r;
        if (s === i) break;
        const tk = this.k[s], tv = this.v[s];
        this.k[s] = this.k[i]; this.v[s] = this.v[i];
        this.k[i] = tk; this.v[i] = tv;
        i = s;
      }
    }
    this.last = { key: rk, val: rv };
    return this.last;
  }
}

/* ----------------------------------------------------------- geometry */
/* 16 grid moves, same order as bonita/router.py OFFSETS ([dc, dr]) */
const OFFSETS = [[1,0],[2,1],[1,1],[1,2],[0,1],[-1,2],[-1,1],[-2,1],
                 [-1,0],[-2,-1],[-1,-1],[-1,-2],[0,-1],[1,-2],[1,-1],[2,-1]];
const NDIR = 16, START_DIR = 16, NSTATE = 17;

/* Segment AB intersects segment CD?  Touching counts (conservative, as in
   geo._segments_intersect -- the boat never gets the benefit of the doubt). */
function segmentsIntersect(ax, ay, bx, by, cx, cy, dx, dy) {
  const o = (px, py, qx, qy, rx, ry) => (qx - px) * (ry - py) - (qy - py) * (rx - px);
  const o1 = o(ax, ay, bx, by, cx, cy), o2 = o(ax, ay, bx, by, dx, dy);
  const o3 = o(cx, cy, dx, dy, ax, ay), o4 = o(cx, cy, dx, dy, bx, by);
  return (o1 * o2) <= 0 && (o3 * o4) <= 0;
}

function pointInPoly(px, py, poly) {
  let inside = false;
  const n = poly.length;
  for (let i = 0, j = n - 1; i < n; j = i++) {
    const yi = poly[i][1], yj = poly[j][1];
    if ((yi > py) !== (yj > py)) {
      const xi = poly[i][0], xj = poly[j][0];
      const x = xi + (py - yi) * (xj - xi) / (yj - yi || 1e-12);
      if (px < x) inside = !inside;
    }
  }
  return inside;
}

class Grid {
  constructor(data, resM, opts) {
    opts = opts || {};
    const P = data.proj, D = data.domain;
    this.data = data;
    this.res = resM;
    this.lat0 = P.lat0; this.lon0 = P.lon0;
    this.mLat = P.mPerDegLat; this.mLon = P.mPerDegLon;

    const x0 = (D.lon_min - this.lon0) * this.mLon;
    const x1 = (D.lon_max - this.lon0) * this.mLon;
    const y0 = (D.lat_min - this.lat0) * this.mLat;
    const y1 = (D.lat_max - this.lat0) * this.mLat;
    this.nx = Math.ceil((x1 - x0) / resM);
    this.ny = Math.ceil((y1 - y0) / resM);
    this.ox = x0; this.oy = y0;
    const N = this.nx * this.ny;
    this.N = N;

    this.X = new Float32Array(N); this.Y = new Float32Array(N);
    this.LON = new Float32Array(N); this.LAT = new Float32Array(N);
    for (let r = 0; r < this.ny; r++) {
      const y = y0 + (r + 0.5) * resM;
      for (let c = 0; c < this.nx; c++) {
        const i = r * this.nx + c, x = x0 + (c + 0.5) * resM;
        this.X[i] = x; this.Y[i] = y;
        this.LON[i] = x / this.mLon + this.lon0;
        this.LAT[i] = y / this.mLat + this.lat0;
      }
    }

    this.land = new Uint8Array(N);
    this.height = new Float32Array(N);
    for (const L of data.land) {
      for (let i = 0; i < N; i++) {
        if (pointInPoly(this.LON[i], this.LAT[i], L.poly)) {
          this.land[i] = 1;
          if (L.height > this.height[i]) this.height[i] = L.height;
        }
      }
    }
    if (opts.avoidShoal && data.shoal) {
      for (let i = 0; i < N; i++)
        if (pointInPoly(this.LON[i], this.LAT[i], data.shoal)) this.land[i] = 1;
    }

    this.distToLand = this._edt();
    const buf = (opts.shoreBufferM === undefined) ? 60 : opts.shoreBufferM;
    this.water = new Uint8Array(N);
    for (let i = 0; i < N; i++)
      this.water[i] = (!this.land[i] && this.distToLand[i] >= buf) ? 1 : 0;

    /*
      SI 10 obstructions -- hard.  Mirrors geo.BayGrid exactly: a water cell
      within `obstructionMarginM` of a line is closed, and every one of the 16
      grid edges whose segment crosses a line is removed from the graph
      (`blocked[i*16+k]`).  Both are needed: the margin alone leaks under a
      2:1 knight move, the edge test alone lets a route graze the line.
      The lines arrive from Python already sealed into the land polygons.
    */
    this.blocked = new Uint8Array(N * NDIR);
    this.obstructed = new Uint8Array(N);
    this.obstructionLines = [];
    const useObs = opts.obstructions !== false && data.obstructions;
    const margin = opts.obstructionMarginM === undefined
      ? (data.obstructionMarginM || 70) : opts.obstructionMarginM;
    if (useObs) {
      for (const ob of data.obstructions) {
        const pts = ob.line.map(p => [(p[0] - this.lon0) * this.mLon, (p[1] - this.lat0) * this.mLat]);
        this.obstructionLines.push({ key: ob.key, pts });
        for (let q = 1; q < pts.length; q++) {
          const ax = pts[q - 1][0], ay = pts[q - 1][1], bx = pts[q][0], by = pts[q][1];
          const vx = bx - ax, vy = by - ay, L2 = Math.max(vx * vx + vy * vy, 1e-9);
          for (let i = 0; i < N; i++) {
            const px = this.X[i], py = this.Y[i];
            let t = ((px - ax) * vx + (py - ay) * vy) / L2; t = t < 0 ? 0 : t > 1 ? 1 : t;
            if (Math.hypot(px - (ax + t * vx), py - (ay + t * vy)) < margin) this.obstructed[i] = 1;
            for (let k = 0; k < NDIR; k++) {
              if (segmentsIntersect(px, py, px + OFFSETS[k][0] * resM, py + OFFSETS[k][1] * resM,
                                    ax, ay, bx, by)) this.blocked[i * NDIR + k] = 1;
            }
          }
        }
      }
      for (let i = 0; i < N; i++) if (this.obstructed[i]) this.water[i] = 0;
    }

    /* shipping zones (soft): index into data.shippingZones, -1 = none */
    this.zone = new Int16Array(N).fill(-1);
    const zones = data.shippingZones || [];
    this.zoneWeight = new Float64Array(zones.length + 1).fill(1);
    for (let k = 0; k < zones.length; k++) {
      this.zoneWeight[k] = zones[k].weight;
      for (let i = 0; i < N; i++)
        if (this.zone[i] < 0 && pointInPoly(this.LON[i], this.LAT[i], zones[k].poly)) this.zone[i] = k;
    }

    this.nWater = 0;
    for (let i = 0; i < N; i++) if (this.water[i]) this.nWater++;
  }

  /* does the straight track through xy points (metres) cross an obstruction line? */
  trackCrossesObstruction(xy) {
    for (let s = 1; s < xy.length; s++) {
      for (const L of this.obstructionLines) {
        for (let q = 1; q < L.pts.length; q++) {
          if (segmentsIntersect(xy[s - 1][0], xy[s - 1][1], xy[s][0], xy[s][1],
                                L.pts[q - 1][0], L.pts[q - 1][1], L.pts[q][0], L.pts[q][1])) return true;
        }
      }
    }
    return false;
  }

  /*
    Exact Euclidean distance to the nearest land cell, in metres
    (Felzenszwalb & Huttenlocher's O(N) squared-distance transform, run down
    the columns then across the rows).

    This has to be exact, not a chamfer approximation: it feeds the shoreline
    current relief, and a 4% error on diagonals moves the current field by
    ~0.3 kt near shore -- which compounds into minutes over a four-hour beat.
    The Python side uses scipy's exact transform, and the two must agree.
  */
  _edt() {
    const nx = this.nx, ny = this.ny, res = this.res, N = this.N;
    const INF = 1e20;
    const f = new Float64Array(Math.max(nx, ny));
    const dsq = new Float64Array(N);
    for (let i = 0; i < N; i++) dsq[i] = this.land[i] ? 0 : INF;

    const v = new Int32Array(Math.max(nx, ny));
    const z = new Float64Array(Math.max(nx, ny) + 1);
    const out = new Float64Array(Math.max(nx, ny));

    const dt1d = (n) => {                     // in-place on f -> out
      let k = 0;
      v[0] = 0; z[0] = -INF; z[1] = INF;
      for (let q = 1; q < n; q++) {
        let s;
        for (;;) {
          s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
          if (s > z[k]) break;
          k--;
        }
        k++; v[k] = q; z[k] = s; z[k + 1] = INF;
      }
      k = 0;
      for (let q = 0; q < n; q++) {
        while (z[k + 1] < q) k++;
        const d = q - v[k];
        out[q] = d * d + f[v[k]];
      }
    };

    for (let c = 0; c < nx; c++) {            // columns
      for (let r = 0; r < ny; r++) f[r] = dsq[r * nx + c];
      dt1d(ny);
      for (let r = 0; r < ny; r++) dsq[r * nx + c] = out[r];
    }
    for (let r = 0; r < ny; r++) {            // rows
      const base = r * nx;
      for (let c = 0; c < nx; c++) f[c] = dsq[base + c];
      dt1d(nx);
      for (let c = 0; c < nx; c++) dsq[base + c] = out[c];
    }
    const d = new Float32Array(N);
    for (let i = 0; i < N; i++) d[i] = Math.sqrt(dsq[i]) * res;
    return d;
  }

  idx(lon, lat) {
    const x = (lon - this.lon0) * this.mLon, y = (lat - this.lat0) * this.mLat;
    let c = Math.round((x - this.ox) / this.res - 0.5);
    let r = Math.round((y - this.oy) / this.res - 0.5);
    c = Math.max(0, Math.min(this.nx - 1, c));
    r = Math.max(0, Math.min(this.ny - 1, r));
    return r * this.nx + c;
  }
  /* containing cell if it is navigable, else the nearest cell that is */
  idxWater(lon, lat) {
    const i = this.idx(lon, lat);
    if (this.water[i]) return i;
    const nx = this.nx, ny = this.ny;
    const r0 = (i / nx) | 0, c0 = i - r0 * nx;
    for (let rad = 1; rad <= 12; rad++) {
      let best = -1, bd = Infinity;
      for (let dr = -rad; dr <= rad; dr++) for (let dc = -rad; dc <= rad; dc++) {
        if (Math.max(Math.abs(dr), Math.abs(dc)) !== rad) continue;
        const r = r0 + dr, c = c0 + dc;
        if (r < 0 || r >= ny || c < 0 || c >= nx) continue;
        const j = r * nx + c;
        if (!this.water[j]) continue;
        const d = dr * dr + dc * dc;
        if (d < bd) { bd = d; best = j; }
      }
      if (best >= 0) return best;
    }
    return -1;
  }
  nearestWater(lon, lat) {
    let best = -1, bd = Infinity;
    const x = (lon - this.lon0) * this.mLon, y = (lat - this.lat0) * this.mLat;
    for (let i = 0; i < this.N; i++) {
      if (!this.water[i]) continue;
      const dx = this.X[i] - x, dy = this.Y[i] - y, d = dx * dx + dy * dy;
      if (d < bd) { bd = d; best = i; }
    }
    return best;
  }
  /*
    Geodesic distance through water from a seed, in metres.

    D must be Float64, not Float32: the heap keys are float64, and storing a
    rounded float32 back into D makes the `d > D[i]` staleness test reject
    valid pops, which silently truncates the search.  That bug reached only
    772 of 3943 water cells and corrupted the current field's interpolation
    weights.
  */
  waterDistance(lon, lat) {
    const N = this.N, nx = this.nx, ny = this.ny, res = this.res;
    const D = new Float64Array(N).fill(Infinity);
    const s = this.nearestWater(lon, lat);
    if (s < 0) return D;
    D[s] = 0;
    const heap = new MinHeap(1024);
    heap.push(0, s);
    const dr = [-1, 1, 0, 0, -1, -1, 1, 1], dc = [0, 0, -1, 1, -1, 1, -1, 1];
    const w = [1, 1, 1, 1, Math.SQRT2, Math.SQRT2, Math.SQRT2, Math.SQRT2];
    let it;
    while ((it = heap.pop())) {
      const i = it.val, d = it.key;
      if (d > D[i]) continue;
      const r = (i / nx) | 0, c = i - r * nx;
      for (let k = 0; k < 8; k++) {
        const rr = r + dr[k], cc = c + dc[k];
        if (rr < 0 || rr >= ny || cc < 0 || cc >= nx) continue;
        const j = rr * nx + cc;
        if (!this.water[j]) continue;
        const nd = d + w[k] * res;
        if (nd < D[j]) { D[j] = nd; heap.push(nd, j); }
      }
    }
    return D;
  }
}

/* ----------------------------------------------------------- currents */
class RefSeries {
  constructor(t, v) { this.t = t; this.v = v; }
  at(h) {
    let x = h % 24; if (x < 0) x += 24;
    const t = this.t, v = this.v, n = t.length;
    if (x <= t[0]) {
      const span = t[0] + 24 - t[n - 1];
      const f = (x + 24 - t[n - 1]) / span;
      return v[n - 1] + f * (v[0] - v[n - 1]);
    }
    if (x >= t[n - 1]) {
      const span = t[0] + 24 - t[n - 1];
      const f = (x - t[n - 1]) / span;
      return v[n - 1] + f * (v[0] - v[n - 1]);
    }
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (t[m] <= x) lo = m; else hi = m; }
    const f = (x - t[lo]) / (t[hi] - t[lo]);
    return v[lo] + f * (v[hi] - v[lo]);
  }
}

class CurrentField {
  /*
    Two tiers.

      measured : the 8 real NOAA current-prediction stations inside the race
                 area.  Each brings its own harmonic series and its own mean
                 flood/ebb set, straight from CO-OPS.
      derived  : water NOAA does not instrument -- Bonita, the Marin shore,
                 Raccoon Strait, the Berkeley Circle, and the two eddies that
                 decide this race.  Each scales a named base station's signed
                 series by a hand amplitude and set.  This tier is the
                 editorial part of the model; it is also the tier to tune.

    Both tiers blend by inverse-distance weighting over distance THROUGH
    WATER, so Sausalito never borrows from the Cityfront across the headlands.
  */
  constructor(grid, data, series, opts) {
    opts = opts || {};
    const p = data.currentParams;
    this.grid = grid;
    this.series = series;                  // {stationId: RefSeries}
    this.strength = opts.strength === undefined ? 1 : opts.strength;

    this.zones = [];
    for (const st of data.stations) {
      if (!series[st.id]) continue;
      this.zones.push({ name: st.name, kind: 'measured', base: st.id, lon: st.lon,
                        lat: st.lat, ampFlood: 1, ampEbb: 1, dirFlood: st.dirFlood,
                        dirEbb: st.dirEbb, lagMin: 0 });
    }
    for (const z of data.derivedZones) {
      if (!series[z.base]) continue;
      this.zones.push(Object.assign({ kind: 'derived' }, z));
    }
    if (!this.zones.length) throw new Error('no current stations available');

    const N = grid.N;
    const W = [], sum = new Float64Array(N);
    for (const z of this.zones) {
      const D = grid.waterDistance(z.lon, z.lat);
      const w = new Float32Array(N);
      for (let i = 0; i < N; i++) {
        w[i] = isFinite(D[i]) ? 1 / Math.pow(D[i] + p.idwD0, p.idwPower) : 0;
        sum[i] += w[i];
      }
      W.push(w);
    }
    for (let k = 0; k < W.length; k++)
      for (let i = 0; i < N; i++) W[k][i] = sum[i] > 0 ? W[k][i] / sum[i] : 0;
    this.w = W;

    this.relief = new Float32Array(N);
    for (let i = 0; i < N; i++)
      this.relief[i] = p.shoreReliefFloor + (1 - p.shoreReliefFloor) *
        Math.tanh(grid.distToLand[i] / p.shoreReliefL);

    this._u = new Float32Array(N);
    this._v = new Float32Array(N);
  }

  nZonesMeasured() { return this.zones.filter(z => z.kind === 'measured').length; }

  at(t) {
    const N = this.grid.N, u = this._u, v = this._v;
    u.fill(0); v.fill(0);
    for (let k = 0; k < this.zones.length; k++) {
      const z = this.zones[k];
      const s = this.series[z.base].at(t - z.lagMin / 60);
      const flood = s >= 0;
      const mag = flood ? z.ampFlood * s : z.ampEbb * (-s);
      const ang = (flood ? z.dirFlood : z.dirEbb) * D2R;
      const uz = mag * Math.sin(ang), vz = mag * Math.cos(ang);
      const w = this.w[k];
      for (let i = 0; i < N; i++) { u[i] += uz * w[i]; v[i] += vz * w[i]; }
    }
    const f = this.strength;
    for (let i = 0; i < N; i++) { const g = this.relief[i] * f; u[i] *= g; v[i] *= g; }
    return { u, v };
  }

  sample(t, lon, lat) {
    const { u, v } = this.at(t), i = this.grid.idxWater(lon, lat);
    return { u: u[i], v: v[i], kt: Math.hypot(u[i], v[i]),
             set: ((Math.atan2(u[i], v[i]) * R2D) % 360 + 360) % 360 };
  }
}

/* --------------------------------------------------------------- wind */
function seaBreeze(p, t) {
  const { baseKt, peakKt, tOnset, tPeak, tHold, tOff } = p;
  if (t < tOnset) return baseKt;
  if (t < tPeak) {
    const f = (t - tOnset) / Math.max(tPeak - tOnset, 1e-6);
    return baseKt + (peakKt - baseKt) * 0.5 * (1 - Math.cos(Math.PI * Math.min(Math.max(f, 0), 1)));
  }
  if (t < tHold) return peakKt;
  if (t < tOff) {
    const f = (t - tHold) / Math.max(tOff - tHold, 1e-6);
    return peakKt + (baseKt - peakKt) * 0.5 * (1 - Math.cos(Math.PI * Math.min(Math.max(f, 0), 1)));
  }
  return baseKt;
}

const DEFAULT_BREEZE = { baseKt: 8, peakKt: 19, tOnset: 10, tPeak: 15.5, tHold: 18, tOff: 21 };

class WindField {
  /*
    Mean strength comes from `profile.mean(t)` -- either the parametric sea
    breeze, or a forecast-driven profile anchored to live anemometer
    observations (see Live.buildProfile).  The spatial structure (Gate
    acceleration, the slot, the ray-marched terrain shadow) multiplies it and
    is normalised to unit mean over water, so `mean(t)` means what it says.
  */
  constructor(grid, data, breeze, opts) {
    opts = opts || {};
    const p = data.windParams;
    this.grid = grid;
    this.breeze = Object.assign({}, DEFAULT_BREEZE, breeze || {});
    this.profile = (opts.profile && typeof opts.profile.mean === 'function')
      ? opts.profile : { mean: (t) => seaBreeze(this.breeze, t) };
    this.dirOffset = 0;          // live veer correction, radians
    this.speedScale = opts.speedScale === undefined ? 1 : opts.speedScale;
    const N = grid.N;
    const dirShift = (opts.dirShiftDeg || 0);

    const pts = data.windPoints;
    const mult = new Float32Array(N), cx = new Float32Array(N), sx = new Float32Array(N);
    const sumw = new Float64Array(N);
    for (const q of pts) {
      const px = (q.lon - grid.lon0) * grid.mLon, py = (q.lat - grid.lat0) * grid.mLat;
      const a = (q.dir + dirShift) * D2R, ca = Math.cos(a), sa = Math.sin(a);
      for (let i = 0; i < N; i++) {
        const d = Math.hypot(grid.X[i] - px, grid.Y[i] - py);
        const w = 1 / Math.pow(d + p.idwD0, p.idwPower);
        sumw[i] += w; mult[i] += w * q.mult; cx[i] += w * ca; sx[i] += w * sa;
      }
    }
    this.mult = mult; this.dirFrom = new Float32Array(N);
    for (let i = 0; i < N; i++) {
      mult[i] /= sumw[i];
      this.dirFrom[i] = Math.atan2(sx[i] / sumw[i], cx[i] / sumw[i]);
    }

    /* directional terrain shadow, ray-marched upwind */
    const lee = new Float32Array(N).fill(1);
    const nSteps = Math.max(1, Math.floor(p.leeMax / grid.res));
    for (let i = 0; i < N; i++) {
      const ux = Math.sin(this.dirFrom[i]), uy = Math.cos(this.dirFrom[i]);
      let hitD = -1, hitH = 0;
      for (let s = 1; s <= nSteps; s++) {
        const d = s * grid.res;
        let c = Math.round((grid.X[i] + ux * d - grid.ox) / grid.res - 0.5);
        let r = Math.round((grid.Y[i] + uy * d - grid.oy) / grid.res - 0.5);
        if (c < 0) c = 0; else if (c >= grid.nx) c = grid.nx - 1;
        if (r < 0) r = 0; else if (r >= grid.ny) r = grid.ny - 1;
        const j = r * grid.nx + c;
        if (grid.land[j]) { hitD = d; hitH = grid.height[j]; break; }
      }
      if (hitD > 0 && hitH > 0) {
        const reach = Math.min(p.leePerHeight * hitH, p.leeMax);
        const deficit = p.leeDeficit * Math.min(Math.max(hitH / p.leeHRef, 0), 1);
        const frac = Math.min(Math.max(hitD / Math.max(reach, 1), 0), 1);
        lee[i] = 1 - deficit * (1 - Math.pow(frac, p.leeExp));
      }
    }
    this.lee = lee;

    let acc = 0, cnt = 0;
    for (let i = 0; i < N; i++) if (grid.water[i]) { acc += mult[i] * lee[i]; cnt++; }
    const mean = acc / Math.max(cnt, 1);
    this.struct = new Float32Array(N);
    for (let i = 0; i < N; i++) this.struct[i] = mult[i] * lee[i] / mean;
    this._s = new Float32Array(N);
  }
  at(t) {
    const s = this.meanAt(t);
    const out = this._s;
    for (let i = 0; i < this.grid.N; i++) out[i] = s * this.struct[i];
    /* dirOffset is a fixed rotation (manual wind); the profile may add a
       time-dependent one (live veer, relaxed -- see Live.buildProfile).
       The returned twd buffer is reused between calls, like tws. */
    const off = this.dirOffset + this.dirShiftAt(t);
    if (off !== 0) {
      if (!this._dirShifted) this._dirShifted = new Float32Array(this.grid.N);
      for (let i = 0; i < this.grid.N; i++)
        this._dirShifted[i] = this.dirFrom[i] + off;
      return { tws: out, twd: this._dirShifted };
    }
    return { tws: out, twd: this.dirFrom };
  }
  meanAt(t) { return this.profile.mean(t) * this.speedScale; }
  dirShiftAt(t) {
    return (this.profile && typeof this.profile.dirShift === 'function')
      ? this.profile.dirShift(t) : 0;
  }
  sample(t, lon, lat) {
    const w = this.at(t), i = this.grid.idxWater(lon, lat);
    return { kt: w.tws[i], from: ((w.twd[i] * R2D) % 360 + 360) % 360 };
  }
  /* structural multiplier at a point -- used to turn a local observation into
     a domain-mean estimate */
  structAt(lon, lat) { return this.struct[this.grid.idxWater(lon, lat)]; }
  /* UNROTATED model direction at a point, radians FROM.  Callers use it to
     compute an offset against the model; including the current offset here
     made a re-applied correction measure against its own previous rotation
     and toggle between full and zero on every GPS fix. */
  dirAt(lon, lat) { return this.dirFrom[this.grid.idxWater(lon, lat)]; }
}

/* -------------------------------------------------------------- polar */
class Polar {
  constructor(data, scale) {
    this.tws = data.polar.tws; this.twa = data.polar.twa;
    this.scale = scale === undefined ? 1 : scale;
    this.sp = data.polar.speed;
    this.name = data.polar.name;
    /* dense lookup table: the router evaluates this ~4000 times per cell, so
       the bilinear walk below is precomputed once onto a regular grid. */
    this.LUT_NA = 181;                       // TWA 0..180 in 1 deg
    this.LUT_NW = 121;                       // TWS 0..30 in 0.25 kt
    this.lut = new Float32Array(this.LUT_NA * this.LUT_NW);
    for (let iw = 0; iw < this.LUT_NW; iw++) {
      const w = iw * 0.25;
      for (let ia = 0; ia < this.LUT_NA; ia++)
        this.lut[iw * this.LUT_NA + ia] = this._interp(ia, w);
    }
  }
  /* fast path used by the router */
  fast(twaDeg, tws) {
    let a = twaDeg % 360; if (a < 0) a += 360; if (a > 180) a = 360 - a;
    let w = tws * 4; if (w < 0) w = 0; else if (w > 120) w = 120;
    const iw = w | 0, fw = w - iw, ia = a | 0, fa = a - ia;
    const iw2 = iw < 120 ? iw + 1 : iw, ia2 = ia < 180 ? ia + 1 : ia;
    const L = this.lut, NA = this.LUT_NA;
    const v = L[iw * NA + ia] * (1 - fw) * (1 - fa) + L[iw2 * NA + ia] * fw * (1 - fa)
            + L[iw * NA + ia2] * (1 - fw) * fa + L[iw2 * NA + ia2] * fw * fa;
    return v * this.scale;
  }
  at(twaDeg, tws) { return this.fast(twaDeg, tws); }
  _interp(twaDeg, tws) {
    let a = Math.abs(((twaDeg + 180) % 360 + 360) % 360 - 180);
    const T = this.tws, A = this.twa;
    let w = Math.min(Math.max(tws, T[0]), T[T.length - 1]);
    let i = 0; while (i < T.length - 2 && T[i + 1] < w) i++;
    let j = 0; while (j < A.length - 2 && A[j + 1] < a) j++;
    const fw = (w - T[i]) / (T[i + 1] - T[i]);
    const fa = (a - A[j]) / (A[j + 1] - A[j]);
    const s = this.sp;
    const v = s[i][j] * (1 - fw) * (1 - fa) + s[i + 1][j] * fw * (1 - fa)
            + s[i][j + 1] * (1 - fw) * fa + s[i + 1][j + 1] * fw * fa;
    return Math.max(v, 0);
  }
  bestVMG(tws, upwind) {
    let bestA = 0, bestV = 0, best = upwind ? -1e9 : 1e9;
    for (let a = 20; a <= 180; a += 0.5) {
      const v = this.at(a, tws), g = v * Math.cos(a * D2R);
      if (upwind ? g > best : g < best) { best = g; bestA = a; bestV = v; }
    }
    return { twa: bestA, bs: bestV, vmg: Math.abs(best) };
  }
}

/*
  Sea state.  Mirrors bonita/polars.py chop_factor exactly.

    wind-against-tide : k is the fractional loss per knot of fully opposed
                        current AT THE BEAT ANGLE.  3 kt dead against in 20 kt
                        costs ~19% beating -- what the Gate does to this boat.
                        A run through the same chop keeps `runFloor` of that,
                        because the bow stops and the kite collapses.
    wind-only chop    : the Bay makes its own slop above ~18 kt regardless of
                        the tide.
*/
function chopFactor(twdFrom, tws, cu, cv, cogDeg, k, opts) {
  k = (k === undefined) ? 0.045 : k;
  opts = opts || {};
  const twaCut = opts.twaCut === undefined ? 95 : opts.twaCut;
  const beatRef = opts.beatRef === undefined ? 40 : opts.beatRef;
  const runFloor = opts.runFloor === undefined ? 0.25 : opts.runFloor;
  const windK = opts.windChopK === undefined ? 0.05 : opts.windChopK;
  const windT = opts.windChopTws === undefined ? 18 : opts.windChopTws;

  const wx = -Math.sin(twdFrom), wy = -Math.cos(twdFrom);
  const cs = Math.hypot(cu, cv);
  let opp = 0;
  if (cs > 1e-6) opp = Math.max(0, Math.min(1, -((cu / cs) * wx + (cv / cs) * wy)));
  const twa = Math.abs(((cogDeg - twdFrom * R2D + 180) % 360 + 360) % 360 - 180);
  let w = Math.max(0, Math.min(1, (twaCut - twa) / Math.max(twaCut - beatRef, 1e-6)));
  if (w < runFloor) w = runFloor;
  const gain = Math.max(0.4, Math.min(1.6, tws / 15));
  let loss = k * gain * w * opp * cs;
  loss += windK * w * Math.max(0, Math.min(1, (tws - windT) / 10));
  return Math.max(0.55, Math.min(1, 1 - loss));
}

/* ------------------------------------------------------------- router */

class Router {
  constructor(grid, windField, currentField, polar, opts) {
    opts = opts || {};
    this.g = grid; this.wf = windField; this.cf = currentField; this.polar = polar;
    this.dtBin = (opts.timeBinMin || 10) / 60;
    this.tackPen = (opts.tackPenaltyS === undefined ? 14 : opts.tackPenaltyS) / 3600;
    this.gybePen = (opts.gybePenaltyS === undefined ? 9 : opts.gybePenaltyS) / 3600;
    this.minSog = opts.minSogKt || 0.15;
    this.chopK = opts.chopK === undefined ? 0.045 : opts.chopK;
    /* An edge is priced with the field frozen at its departure bin.  Allowing
       an edge to run for many bins lets the router exploit stale water -- it
       will "sail" one cell at 0.18 kt for 100 minutes through a flood that has
       long since eased.  That is waiting, not sailing, and this formulation
       cannot represent waiting, so it is refused instead.  See holdAdvice(). */
    this.maxEdgeH = (opts.maxEdgeBins === undefined ? 3 : opts.maxEdgeBins) * this.dtBin;
    this.currentScale = opts.currentScale === undefined ? 1 : opts.currentScale;
    /* Soft shipping-lane cost (SI 9.1): time spent on an edge that enters a
       zone cell is multiplied by the zone weight in the search objective
       only; real time is carried alongside and is what gets reported.
       laneCost:false keeps the crossing flags and charges nothing. */
    this.laneCost = opts.laneCost !== false;
    this.zoneW = new Float64Array(grid.zoneWeight.length).fill(1);
    if (this.laneCost) this.zoneW.set(grid.zoneWeight);

    const res = grid.res;
    this.legNm = new Float64Array(NDIR);
    this.gx = new Float64Array(NDIR); this.gy = new Float64Array(NDIR);
    this.cog = new Float64Array(NDIR);
    this.stepIdx = new Int32Array(NDIR);
    for (let k = 0; k < NDIR; k++) {
      const dx = OFFSETS[k][0] * res, dy = OFFSETS[k][1] * res;
      const L = Math.hypot(dx, dy);
      this.legNm[k] = L / 1852;
      this.gx[k] = dx / L; this.gy[k] = dy / L;
      this.cog[k] = ((Math.atan2(dx, dy) * R2D) % 360 + 360) % 360;
      this.stepIdx[k] = OFFSETS[k][1] * grid.nx + OFFSETS[k][0];
    }
    const NS = opts.headingSamples || 180;
    this.th = new Float64Array(NS); this.sinT = new Float64Array(NS);
    this.cosT = new Float64Array(NS); this.thDeg = new Float64Array(NS);
    for (let i = 0; i < NS; i++) {
      const a = 2 * Math.PI * i / NS;
      this.th[i] = a; this.sinT[i] = Math.sin(a); this.cosT[i] = Math.cos(a);
      this.thDeg[i] = a * R2D;
    }
    this.NS = NS;
    this._vx = new Float64Array(NS); this._vy = new Float64Array(NS);
    this.fieldCache = new Map();
    this.edgeCache = new Map();
    this._edgeBuf = new Float64Array(NDIR * 2);
  }

  fields(t) {
    const b = Math.floor(t / this.dtBin);
    let hit = this.fieldCache.get(b);
    if (!hit) {
      const tb = (b + 0.5) * this.dtBin;
      const w = this.wf.at(tb), c = this.cf.at(tb);
      hit = { bin: b, tws: Float32Array.from(w.tws), twd: Float32Array.from(w.twd),
              cu: Float32Array.from(c.u), cv: Float32Array.from(c.v) };
      if (this.currentScale !== 1) {
        for (let i = 0; i < hit.cu.length; i++) { hit.cu[i] *= this.currentScale; hit.cv[i] *= this.currentScale; }
      }
      this.fieldCache.set(b, hit);
    }
    return hit;
  }

  /* SOG on each of the 16 grid directions from cell i at time t */
  edgeSpeeds(i, t) {
    const F = this.fields(t);
    const key = i * 100000 + (F.bin + 50000);
    let hit = this.edgeCache.get(key);
    if (hit) return hit;

    const w = F.tws[i], wd = F.twd[i], cu = F.cu[i], cv = F.cv[i];
    const NS = this.NS, vx = this._vx, vy = this._vy;
    const wdDeg = wd * R2D;
    for (let s = 0; s < NS; s++) {
      const twa = this.thDeg[s] - wdDeg;
      let bs = this.polar.at(twa, w);
      bs *= chopFactor(wd, w, cu, cv, this.thDeg[s], this.chopK);
      vx[s] = bs * this.sinT[s] + cu;
      vy[s] = bs * this.cosT[s] + cv;
    }
    const sog = new Float64Array(NDIR);
    const hd = new Float64Array(NDIR);
    for (let k = 0; k < NDIR; k++) {
      const gx = this.gx[k], gy = this.gy[k];
      let best = -Infinity, bestH = NaN;
      let prevCross = vx[NS - 1] * gy - vy[NS - 1] * gx;
      let prevDot = vx[NS - 1] * gx + vy[NS - 1] * gy;
      for (let s = 0; s < NS; s++) {
        const cross = vx[s] * gy - vy[s] * gx;
        const dot = vx[s] * gx + vy[s] * gy;
        if ((prevCross <= 0 && cross > 0) || (prevCross >= 0 && cross < 0)) {
          const den = prevCross - cross;
          let f = Math.abs(den) > 1e-12 ? prevCross / den : 0;
          f = f < 0 ? 0 : f > 1 ? 1 : f;
          const v = prevDot * (1 - f) + dot * f;
          if (v > best) {
            best = v;
            bestH = ((this.thDeg[(s - 1 + NS) % NS] + f * (360 / NS)) % 360 + 360) % 360;
          }
        }
        prevCross = cross; prevDot = dot;
      }
      /* refuse edges that would take longer than the field is valid for */
      if (best > this.minSog && this.legNm[k] / best > this.maxEdgeH) best = -1;
      sog[k] = best > this.minSog ? best : -1;
      hd[k] = bestH;
    }
    hit = { sog, heading: hd, twd: wd, tws: w, cu, cv };
    this.edgeCache.set(key, hit);
    return hit;
  }

  turnPenalty(iIn, iOut, wd) {
    if (iIn === START_DIR || iIn === iOut) return 0;
    const wdDeg = wd * R2D;
    const ai = ((this.cog[iIn] - wdDeg + 180) % 360 + 360) % 360 - 180;
    const ao = ((this.cog[iOut] - wdDeg + 180) % 360 + 360) % 360 - 180;
    if (ai * ao < 0) {
      /* which axis it crossed is what tells a tack from a gybe: the two angles
         from the wind sum to < 180 through head-to-wind, > 180 through dead
         downwind. */
      return (Math.abs(ai) + Math.abs(ao) < 180) ? this.tackPen : this.gybePen;
    }
    return 0;
  }

  route(startLL, endLL, t0, opts) {
    opts = opts || {};
    const g = this.g, nx = g.nx, ny = g.ny, N = g.N;
    const forbid = opts.forbid || null;
    const okCell = (i) => g.water[i] && (!forbid || !forbid[i]);
    const s = g.nearestWater(startLL[0], startLL[1]);
    const e = g.nearestWater(endLL[0], endLL[1]);
    if (s < 0 || e < 0) return { ok: false, note: 'endpoint not navigable' };
    /* Do NOT relocate an endpoint into the corridor.  Silently moving the
       start to the nearest in-lane cell made constrained lanes come out FASTER
       than the unconstrained optimum, which is impossible, and would have sent
       the boat to the wrong side of the Bay. */
    if (forbid && !okCell(s))
      return { ok: false, note: 'you are outside this lane' };
    if (forbid && !okCell(e))
      return { ok: false, note: 'the mark is outside this lane' };

    const best = new Float64Array(N * NSTATE).fill(Infinity);   // cost labels
    const treal = new Float64Array(N * NSTATE).fill(Infinity);  // real time at label
    const prev = new Int32Array(N * NSTATE).fill(-1);
    best[s * NSTATE + START_DIR] = t0;
    treal[s * NSTATE + START_DIR] = t0;
    const heap = new MinHeap(1 << 15);
    heap.push(t0, s * NSTATE + START_DIR);
    const blocked = g.blocked, zone = g.zone, zw = this.zoneW, nz = zw.length - 1;
    let goal = -1, goalT = Infinity, goalK = Infinity;
    let pops = 0;
    const cap = opts.maxPops || 4e6;

    let it;
    while ((it = heap.pop())) {
      const st = it.val, k = it.key;
      if (k > best[st] + 1e-12) continue;
      const cell = (st / NSTATE) | 0, di = st - cell * NSTATE;
      const t = treal[st];
      if (cell === e) { goal = st; goalT = t; goalK = k; break; }
      if (++pops > cap) break;
      const E = this.edgeSpeeds(cell, t);
      const r = (cell / nx) | 0, c = cell - r * nx;
      const bb = cell * NDIR;
      for (let d = 0; d < NDIR; d++) {
        if (blocked[bb + d]) continue;
        const sog = E.sog[d];
        if (sog <= 0) continue;
        const rr = r + OFFSETS[d][1], cc = c + OFFSETS[d][0];
        if (rr < 0 || rr >= ny || cc < 0 || cc >= nx) continue;
        const j = rr * nx + cc;
        if (!okCell(j)) continue;
        const dt = this.legNm[d] / sog + this.turnPenalty(di, d, E.twd);
        const zj = zone[j];
        const nk = k + dt * zw[zj < 0 ? nz : zj];
        const ns = j * NSTATE + d;
        if (nk < best[ns] - 1e-12) { best[ns] = nk; treal[ns] = t + dt; prev[ns] = st; heap.push(nk, ns); }
      }
    }
    if (goal < 0) return { ok: false, note: 'no route found' };

    const path = [];
    let st = goal;
    while (st >= 0) {
      const cell = (st / NSTATE) | 0;
      path.push({ i: cell, t: treal[st], lon: g.LON[cell], lat: g.LAT[cell],
                  dir: st - cell * NSTATE });
      st = prev[st];
    }
    path.reverse();
    return { ok: true, tStart: t0, tEnd: goalT, elapsed: goalT - t0, path,
             costEnd: goalK, lanePenalty: goalK - goalT,
             zoneCrossings: zoneCrossings(g, path) };
  }

  _nearestOk(ll, okCell) {
    const g = this.g;
    const x = (ll[0] - g.lon0) * g.mLon, y = (ll[1] - g.lat0) * g.mLat;
    let best = -1, bd = Infinity;
    for (let i = 0; i < g.N; i++) {
      if (!okCell(i)) continue;
      const dx = g.X[i] - x, dy = g.Y[i] - y, d = dx * dx + dy * dy;
      if (d < bd) { bd = d; best = i; }
    }
    return best;
  }
}

/* Where and when a path is inside a shipping zone (merged runs, path order).
   Mirrors router.zone_crossings. */
function zoneCrossings(grid, path) {
  const zones = grid.data.shippingZones || [];
  const out = [];
  let cur = null;
  for (const p of path) {
    const k = grid.zone[p.i];
    if (cur && k === cur.k) { cur.tOut = p.t; continue; }
    if (cur) { out.push(cur); cur = null; }
    if (k >= 0) cur = { k, key: zones[k].key, zone: zones[k].name, kind: zones[k].kind, tIn: p.t, tOut: p.t };
  }
  if (cur) out.push(cur);
  return out;
}

/* ---------------------------------------------------------- live data */
const COOPS = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter';
const OPENMETEO = 'https://api.open-meteo.com/v1/forecast';

function angDiffDeg(a, b) { return ((a - b + 540) % 360) - 180; }

const Live = {
  cacheKey: (k) => 'bonita.cache.' + k,

  cacheGet(k, maxAgeMs) {
    try {
      const raw = localStorage.getItem(Live.cacheKey(k));
      if (!raw) return null;
      const o = JSON.parse(raw);
      if (maxAgeMs && Date.now() - o.at > maxAgeMs) return null;
      return o;
    } catch (e) { return null; }
  },
  cacheSet(k, value) {
    try { localStorage.setItem(Live.cacheKey(k), JSON.stringify({ at: Date.now(), value })); }
    catch (e) { /* quota or private mode -- caching is best effort */ }
  },

  /* --- NOAA harmonic current predictions, one call per station ---------- */
  async fetchCurrents(dateStr, stationIds, opts) {
    opts = opts || {};
    const ymd = dateStr.replace(/-/g, '');
    const key = 'cur.' + ymd;
    const cached = Live.cacheGet(key, opts.maxAgeMs || 30 * 24 * 3600e3);
    const have = cached && cached.value || {};
    const out = {}, fresh = {};
    let fetched = 0, failed = 0;
    await Promise.all(stationIds.map(async (id) => {
      if (have[id] && !opts.force) { out[id] = have[id]; return; }
      try {
        const url = `${COOPS}?product=currents_predictions&application=bonita_router`
          + `&begin_date=${ymd}&end_date=${ymd}&station=${id}`
          + `&time_zone=lst_ldt&interval=30&units=english&format=json`;
        const r = await fetch(url);
        const j = await r.json();
        const cp = j.current_predictions && j.current_predictions.cp;
        if (!cp || !cp.length) throw new Error('no data');
        const t = [], v = [];
        for (const row of cp) {
          const hm = row.Time.slice(11);
          t.push(+hm.slice(0, 2) + (+hm.slice(3, 5)) / 60);
          v.push(+row.Velocity_Major);
        }
        out[id] = { t, v, dirFlood: +cp[0].meanFloodDir, dirEbb: +cp[0].meanEbbDir };
        fresh[id] = out[id]; fetched++;
      } catch (e) { failed++; if (have[id]) out[id] = have[id]; }
    }));
    if (fetched) Live.cacheSet(key, Object.assign({}, have, fresh));
    return { series: out, fetched, failed, cachedOnly: fetched === 0 };
  },

  /* --- live observed wind from CO-OPS anemometers ----------------------- */
  async fetchWindObs(metStations) {
    const obs = [];
    await Promise.all(metStations.map(async (m) => {
      try {
        const url = `${COOPS}?product=wind&application=bonita_router&date=latest`
          + `&station=${m.id}&time_zone=lst_ldt&units=english&format=json`;
        const r = await fetch(url);
        const j = await r.json();
        const d = j.data && j.data[j.data.length - 1];
        if (!d) throw new Error('no data');
        /* CO-OPS `units=english` is documented as "fahrenheit, feet, KNOTS".
           An earlier version multiplied by 0.868976 (the mph->kt factor) and
           read every anemometer 13% low, which cost ~18 min of elapsed time. */
        const kt = (d.s === '' || d.s == null) ? NaN : +d.s;
        const gust = (d.g === '' || d.g == null) ? NaN : +d.g;
        if (!isFinite(kt)) throw new Error('empty reading');
        obs.push({ id: m.id, name: m.name, lat: m.lat, lon: m.lon,
                   kt, gustKt: gust, from: +d.d, time: d.t });
      } catch (e) { /* station offline; the others still anchor the field */ }
    }));
    if (obs.length) Live.cacheSet('windobs', obs);
    else { const c = Live.cacheGet('windobs', 6 * 3600e3); if (c) return { obs: c.value, stale: true }; }
    return { obs, stale: false };
  },

  /* --- gridded forecast ------------------------------------------------ */
  async fetchWindForecast(dateStr, pts) {
    const key = 'fc.' + dateStr;
    try {
      const lat = pts.map(p => p[1].toFixed(3)).join(',');
      const lon = pts.map(p => p[0].toFixed(3)).join(',');
      const url = `${OPENMETEO}?latitude=${lat}&longitude=${lon}`
        + `&hourly=wind_speed_10m,wind_direction_10m&wind_speed_unit=kn`
        + `&timezone=America%2FLos_Angeles&start_date=${dateStr}&end_date=${dateStr}`;
      const r = await fetch(url);
      let j = await r.json();
      if (!Array.isArray(j)) j = [j];
      const times = j[0].hourly.time.map(s => +s.slice(11, 13) + (+s.slice(14, 16)) / 60);
      const val = { times, speed: j.map(d => d.hourly.wind_speed_10m),
                    dir: j.map(d => d.hourly.wind_direction_10m),
                    lat: pts.map(p => p[1]), lon: pts.map(p => p[0]) };
      Live.cacheSet(key, val);
      return { fc: val, stale: false };
    } catch (e) {
      const c = Live.cacheGet(key, 7 * 24 * 3600e3);
      return c ? { fc: c.value, stale: true } : { fc: null, stale: true };
    }
  },

  /*
    Build the mean-strength profile the wind field runs on.

      forecast  sets the SHAPE through the day (spatial mean of the grid
                points, divided by nothing -- the field structure carries the
                spatial part).
      obs       sets the LEVEL right now: each anemometer reading is divided
                by the field's structural multiplier at that spot to get an
                implied domain mean, and those are averaged.
      the anchor correction relaxes back to the raw forecast over
      `relaxHours`, because an observation tells you about now, not 16:00.
  */
  buildProfile(windField, opts) {
    opts = opts || {};
    const fc = opts.forecast, obsIn = opts.obs || [], tNow = opts.tNow;
    const relax = opts.relaxHours || 3;
    const dom = opts.domain;
    const fallback = (t) => seaBreeze(windField.breeze, t);

    let fcMean = null;
    if (fc && fc.times && fc.times.length) {
      const n = fc.times.length, m = new Float64Array(n);
      for (let k = 0; k < n; k++) {
        let s = 0, c = 0;
        for (let p = 0; p < fc.speed.length; p++) {
          const v = fc.speed[p][k];
          if (v != null && isFinite(v)) { s += v; c++; }
        }
        m[k] = c ? s / c : NaN;
      }
      fcMean = (t) => {
        const T = fc.times;
        if (t <= T[0]) return m[0];
        if (t >= T[n - 1]) return m[n - 1];
        let lo = 0, hi = n - 1;
        while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (T[mid] <= t) lo = mid; else hi = mid; }
        const f = (t - T[lo]) / (T[hi] - T[lo]);
        return m[lo] + f * (m[hi] - m[lo]);
      };
    }
    const base = fcMean || fallback;

    /*
      Only anchor on anemometers that actually sit inside the model domain.
      A station outside it resolves to a clamped edge cell, so its structural
      multiplier is an artefact of the boundary rather than of its real
      shelter -- and dividing a sheltered harbour reading by ~1.0 then pulls
      the whole Bay's wind down by a third.
    */
    const obs = obsIn.filter(o => !dom ||
      (o.lon > dom.lon_min && o.lon < dom.lon_max &&
       o.lat > dom.lat_min && o.lat < dom.lat_max));

    /* the observations describe the moment they were taken, not "now" */
    let tObs = tNow;
    const stamps = obs.map(o => o.time).filter(Boolean);
    if (stamps.length) {
      const hm = stamps[0].slice(11);
      const h = +hm.slice(0, 2) + (+hm.slice(3, 5)) / 60;
      if (isFinite(h)) tObs = h;
    }

    let anchor = 1, impliedNow = null, dirErr = null;
    if (obs.length && tObs != null) {
      const implied = [], dErr = [];
      for (const o of obs) {
        const st = windField.structAt(o.lon, o.lat);
        if (isFinite(st) && st > 0.05 && isFinite(o.kt)) implied.push(o.kt / st);
        const md = windField.dirAt(o.lon, o.lat);
        if (isFinite(o.from) && isFinite(md) && o.kt > 3)
          dErr.push(angDiffDeg(o.from, md * R2D));
      }
      if (implied.length) {
        implied.sort((a, b) => a - b);
        const n = implied.length;
        impliedNow = n % 2 ? implied[(n - 1) / 2]
                           : 0.5 * (implied[n / 2 - 1] + implied[n / 2]);
        const b = base(tObs);
        if (b > 0.5) anchor = impliedNow / b;
      }
      if (dErr.length) {
        let sx = 0, cx = 0;
        for (const d of dErr) { sx += Math.sin(d * D2R); cx += Math.cos(d * D2R); }
        const e = Math.atan2(sx / dErr.length, cx / dErr.length) * R2D;
        if (isFinite(e)) dirErr = e;
      }
    }

    const mean = (t) => {
      const b = base(t);
      if (anchor === 1 || tObs == null) return b;
      return b * (1 + (anchor - 1) * Math.exp(-Math.abs(t - tObs) / relax));
    };
    /*
      The direction error relaxes exactly like the strength anchor.  Applied
      as a fixed rotation of the whole day it turned an 08:24 south-westerly
      land breeze (-38 deg vs the model) into a -38 deg rotation of the 13:00
      sea breeze, and the beat to Bonita came out as four straight boards
      sailed within 8 deg of the real wind.  Radians, added to dirFrom.
    */
    const dirShift = (t) => (dirErr == null || tObs == null) ? 0
      : dirErr * D2R * Math.exp(-Math.abs(t - tObs) / relax);

    return { mean, dirShift, source: fcMean ? 'forecast' : 'parametric',
             anchored: anchor !== 1, anchor, impliedNow, dirErr,
             nObs: obs.length, nDropped: obsIn.length - obs.length, tObs };
  },

  /*
    GPS check: what does the model think you should be doing right now, and
    what are you actually doing?  Returns predicted SOG for your observed COG
    at your observed position, and the ratio.  A persistent ratio away from 1
    means the wind or the polar is off -- the UI offers to scale the wind by
    it, which is the honest interpretation most of the time.
  */
  gpsCheck(router, lon, lat, tNow, cogDeg, sogKt) {
    const g = router.g;
    const i = g.idxWater(lon, lat);
    if (i < 0) return { ok: false, reason: 'off the model grid' };
    const F = router.fields(tNow);
    const w = F.tws[i], wd = F.twd[i], cu = F.cu[i], cv = F.cv[i];
    const NS = router.NS;
    let best = -Infinity;
    const gx = Math.sin(cogDeg * D2R), gy = Math.cos(cogDeg * D2R);
    let prevCross = null, prevDot = null;
    for (let s = 0; s <= NS; s++) {
      const k = s % NS;
      const twa = router.thDeg[k] - wd * R2D;
      let bs = router.polar.at(twa, w);
      bs *= chopFactor(wd, w, cu, cv, router.thDeg[k], router.chopK);
      const vx = bs * router.sinT[k] + cu, vy = bs * router.cosT[k] + cv;
      const cross = vx * gy - vy * gx, dot = vx * gx + vy * gy;
      if (prevCross !== null &&
          ((prevCross <= 0 && cross > 0) || (prevCross >= 0 && cross < 0))) {
        const den = prevCross - cross;
        let f = Math.abs(den) > 1e-12 ? prevCross / den : 0;
        f = f < 0 ? 0 : f > 1 ? 1 : f;
        const v = prevDot * (1 - f) + dot * f;
        if (v > best) best = v;
      }
      prevCross = cross; prevDot = dot;
    }
    if (!isFinite(best) || best <= 0) return { ok: false, reason: 'course not sailable in the model' };
    return { ok: true, predictedSog: best, observedSog: sogKt,
             ratio: sogKt / best, tws: w, twdDeg: ((wd * R2D) % 360 + 360) % 360,
             curKt: Math.hypot(cu, cv),
             curSet: ((Math.atan2(cu, cv) * R2D) % 360 + 360) % 360 };
  },
};

/* ------------------------------------------------------------ helpers */
function trackNm(path, grid) {
  let s = 0;
  for (let k = 1; k < path.length; k++) {
    const a = path[k - 1], b = path[k];
    const dx = (b.lon - a.lon) * grid.mLon, dy = (b.lat - a.lat) * grid.mLat;
    s += Math.hypot(dx, dy);
  }
  return s / 1852;
}

function hhmm(t) {
  let x = ((t % 24) + 24) % 24;
  let h = Math.floor(x), m = Math.round((x - h) * 60);
  if (m === 60) { m = 0; h = (h + 1) % 24; }
  return String(h).padStart(2, '0') + ':' + String(m).padStart(2, '0');
}

function bearing(a, b, grid) {
  const dx = (b[0] - a[0]) * grid.mLon, dy = (b[1] - a[1]) * grid.mLat;
  return ((Math.atan2(dx, dy) * R2D) % 360 + 360) % 360;
}
function distNm(a, b, grid) {
  const dx = (b[0] - a[0]) * grid.mLon, dy = (b[1] - a[1]) * grid.mLat;
  return Math.hypot(dx, dy) / 1852;
}

root.BonitaModel = { Grid, CurrentField, WindField, Polar, Router, RefSeries, Live,
                     MinHeap, seaBreeze, chopFactor, DEFAULT_BREEZE, OFFSETS,
                     segmentsIntersect, zoneCrossings, pointInPoly,
                     trackNm, hhmm, bearing, distNm, angDiffDeg, D2R, R2D };
})(typeof globalThis !== 'undefined' ? globalThis : this);
