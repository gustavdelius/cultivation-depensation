// Check the JavaScript engine against trajectories exported from mizer.
// Run from the project root:  node interactive/validate.js
"use strict";
const path = require("path");
const MODEL = require(path.join(__dirname, "model_data.js"));
const ref = require(path.join(__dirname, "validation", "reference.json"));
const { createEngine } = require(path.join(__dirname, "engine.js"));

let failures = 0;
function check(label, got, want, tol) {
  const err = Math.abs(got - want) / Math.max(Math.abs(want), 1e-12);
  const ok = err <= tol;
  if (!ok) failures++;
  console.log(`${ok ? "ok  " : "FAIL"} ${label}: js=${got.toPrecision(6)} mizer=${want.toPrecision(6)} rel.err=${err.toExponential(1)}`);
}
function maxRelErr(a, b) {
  let m = 0, scale = 0;
  for (let i = 0; i < b.length; i++) scale = Math.max(scale, Math.abs(b[i]));
  for (let i = 0; i < b.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]) / scale);
  return m;
}

const eng = createEngine(MODEL);

// 1. Rates at State RFP
eng.loadState("State_RFP");
const rRFP = eng.rates(eng.state.n, eng.state.n_pp, [0, 0]);
for (const [name, want] of Object.entries(ref.rates_RFP)) {
  const flatWant = [].concat(...(Array.isArray(want[0]) ? want : [want]));
  const err = maxRelErr(rRFP[name], flatWant);
  const ok = err < 1e-9;
  if (!ok) failures++;
  console.log(`${ok ? "ok  " : "FAIL"} rate ${name} at State RFP: max err/scale = ${err.toExponential(2)}`);
}

// 2. R0 at the four unfished states
for (const [st, want] of Object.entries(ref.R0)) {
  eng.loadState(st);
  eng.state.effort.set([0, 0]);
  const d = eng.diagnostics();
  check(`R0 Forage   ${st}`, d.cohort[0].R0, want.Forage, 1e-8);
  check(`R0 Predator ${st}`, d.cohort[1].R0, want.Predator, 1e-8);
}

// 3. Trajectories
for (const key of ["collapse", "balanced", "invade_RP", "rescue_RF"]) {
  const run = ref[key];
  eng.loadState(run.start);
  for (const sd of [].concat(run.seed || [])) {
    for (let j = 0; j < eng.dims.NW; j++) eng.state.n[sd.species * eng.dims.NW + j] = 0;
    eng.seed(sd.species, sd.factor);
  }
  eng.state.effort.set(run.effort);
  const stepsPerYear = Math.round(1 / ref.dt);
  let worst = 0;
  for (let y = 0; y < run.t.length; y++) {
    if (y > 0) for (let k = 0; k < stepsPerYear; k++) eng.step(ref.dt);
    for (const [s, sp] of [[0, "forage"], [1, "predator"]]) {
      const want = run[sp][y];
      const scale = Math.max(...run[sp]);
      worst = Math.max(worst, Math.abs(eng.biomass(s) - want) / scale);
    }
  }
  const last = run.t.length - 1;
  const ok = worst < 1e-6;
  if (!ok) failures++;
  console.log(`${ok ? "ok  " : "FAIL"} trajectory ${key}: max |err|/max = ${worst.toExponential(2)}; ` +
    `final B_forage js=${eng.biomass(0).toFixed(2)} mizer=${run.forage[last].toFixed(2)}, ` +
    `B_pred js=${eng.biomass(1).toFixed(2)} mizer=${run.predator[last].toFixed(2)}`);
}

// 4. Timing
const t0 = Date.now();
eng.loadState("State_RFP");
for (let i = 0; i < 5000; i++) eng.step(0.2);
console.log(`speed: ${(5000 * 0.2 / ((Date.now() - t0) / 1000)).toFixed(0)} model years per second`);

console.log(failures ? `\n${failures} check(s) FAILED` : "\nAll checks passed");
process.exit(failures ? 1 : 0);
