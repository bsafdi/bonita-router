/*
  Route both SI 12 courses with the browser model, using the Python router's
  settings (5-minute bins, 240 headings), and print the numbers the Python
  test compares against.  Run by tests/test_solver.py when node is present:
      node tests/js_parity.js web/model_data.json
*/
const path = require('path');
require(path.join(__dirname, '..', 'web', 'model.js'));
const M = globalThis.BonitaModel;
const DATA = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));
const g = new M.Grid(DATA, 250);
let nb = 0, nobs = 0;
for (let i = 0; i < g.N; i++) {
  if (g.obstructed[i]) nobs++;
  for (let k = 0; k < 16; k++) if (g.blocked[i * 16 + k]) { nb++; break; }
}
const series = {};
for (const [id, s] of Object.entries(DATA.fallbackSeries)) series[id] = new M.RefSeries(s.t, s.v);
const cf = new M.CurrentField(g, DATA, series, {});
const wf = new M.WindField(g, DATA, {});
const R = new M.Router(g, wf, cf, new M.Polar(DATA, 1), { timeBinMin: 5, headingSamples: 240 });
const out = { nWater: g.nWater, nObstructed: nobs, nBlockedCells: nb, courses: {} };
for (const c of Object.keys(DATA.courses)) {
  const mark = DATA.marks[DATA.courses[c].mark];
  const a = R.route(DATA.marks.start, mark, DATA.race.startH);
  const b = R.route(mark, DATA.marks.finish, a.tEnd);
  const xy = a.path.concat(b.path).map(p => [(p.lon - g.lon0) * g.mLon, (p.lat - g.lat0) * g.mLat]);
  out.courses[c] = { tRound: a.tEnd, tFinish: b.tEnd, nOut: a.path.length, nBack: b.path.length,
                     crosses: g.trackCrossesObstruction(xy),
                     zonesOut: a.zoneCrossings.map(z => z.key), zonesBack: b.zoneCrossings.map(z => z.key) };
}
console.log(JSON.stringify(out));
