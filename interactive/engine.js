// Size-spectrum engine: a JavaScript port of mizer's default dynamics for the
// cultivation-depensation model (see export_model.R for the data it consumes).
//
// It reproduces mizer 3.4's `project()` with method = "euler":
//   * encounter and predation rate as the same circular convolutions that
//     mizer evaluates by FFT (wrap-around included),
//   * Holling type II feeding, energy split by the maturity ogive psi,
//   * Beverton-Holt recruitment,
//   * semichemostat resource updated with its exact exponential solution,
//   * implicit first-order upwind transport in size.
// validate.js checks it against trajectories exported from mizer.

(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.SizeSpectrumEngine = factory();
})(this, function () {
  "use strict";

  function createEngine(M) {
    const S = M.species.length;          // number of species
    const NW = M.w.length;               // fish size bins
    const NF = M.w_full.length;          // full size grid (resource + fish)
    const off = M.offset0;               // w[j] == w_full[j + off]

    const flat = (rows) => Float64Array.from([].concat(...rows));
    const sv = flat(M.search_vol), imax = flat(M.intake_max), metab = flat(M.metab);
    const psi = flat(M.psi), mu_b = flat(M.mu_b), Fpe = flat(M.F_per_effort);
    const phiE = flat(M.phi_e), phiP = flat(M.phi_p), ftMask = flat(M.ft_mask);
    const w = Float64Array.from(M.w), dw = Float64Array.from(M.dw);
    const wdwFull = Float64Array.from(M.w_full.map((x, k) => x * M.dw_full[k]));
    const rr = Float64Array.from(M.rr_pp), cc = Float64Array.from(M.cc_pp);
    const theta = M.interaction, thetaR = M.interaction_resource;
    const wMinIdx = M.w_min_idx;
    const wTop = M.w_max_idx.map((i) => Math.min(i + 1, NW)); // support top (count)

    // Work arrays reused every step
    const prey = new Float64Array(NF);
    const Q = new Float64Array(NF);
    const r = {
      encounter: new Float64Array(S * NW), feeding_level: new Float64Array(S * NW),
      e_repro: new Float64Array(S * NW), e_growth: new Float64Array(S * NW),
      pred_rate: new Float64Array(S * NF), pred_mort: new Float64Array(S * NW),
      f_mort: new Float64Array(S * NW), mort: new Float64Array(S * NW),
      rdi: new Float64Array(S), rdd: new Float64Array(S),
      resource_mort: new Float64Array(NF),
    };

    const state = {
      n: new Float64Array(S * NW),
      n_pp: new Float64Array(NF),
      effort: new Float64Array(S),
      t: 0,
    };

    function rates(n, n_pp, effort) {
      for (let i = 0; i < S; i++) {
        // Prey biomass density seen by species i, times w dw (mizerEncounter)
        for (let k = 0; k < NF; k++) prey[k] = thetaR[i] * n_pp[k];
        for (let s = 0; s < S; s++) {
          const th = theta[i][s];
          if (th === 0) continue;
          for (let j = 0; j < NW; j++) prey[j + off] += th * n[s * NW + j];
        }
        for (let k = 0; k < NF; k++) prey[k] *= wdwFull[k];
        // Circular convolution, only at the fish sizes we keep
        const pe = i * NF;
        for (let j = 0; j < NW; j++) {
          const J = j + off;
          let sum = 0;
          for (let k = 0; k < NF; k++) {
            let m = J - k;
            if (m < 0) m += NF;
            sum += phiE[pe + m] * prey[k];
          }
          if (sum < 1e-18) sum = 0;
          const ij = i * NW + j;
          const enc = sv[ij] * sum;
          const f = enc / (enc + imax[ij]);
          r.encounter[ij] = enc;
          r.feeding_level[ij] = f;
          const e = M.alpha[i] * (1 - f) * enc - metab[ij];
          const er = psi[ij] * (e > 0 ? e : 0);
          r.e_repro[ij] = er;
          const g = e - er;
          r.e_growth[ij] = g > 0 ? g : 0;
        }
        // Predation rate exerted by species i on prey of each size
        Q.fill(0);
        for (let j = 0; j < NW; j++) {
          const ij = i * NW + j;
          Q[j + off] = (1 - r.feeding_level[ij]) * sv[ij] * n[ij] * dw[j];
        }
        for (let k = 0; k < NF; k++) {
          let sum = 0;
          for (let J = off; J < NF; J++) {
            let m = k - J;
            if (m < 0) m += NF;
            sum += phiP[pe + m] * Q[J];
          }
          if (sum < 1e-18) sum = 0;
          r.pred_rate[pe + k] = sum * ftMask[pe + k];
        }
      }
      // Mortality and reproduction
      for (let s = 0; s < S; s++) {
        let repro = 0;
        for (let j = 0; j < NW; j++) {
          const sj = s * NW + j;
          let pm = 0;
          for (let i = 0; i < S; i++) pm += theta[i][s] * r.pred_rate[i * NF + j + off];
          r.pred_mort[sj] = pm;
          r.f_mort[sj] = effort[s] * Fpe[sj];
          r.mort[sj] = pm + mu_b[sj] + r.f_mort[sj];
          repro += r.e_repro[sj] * n[sj] * dw[j];
        }
        const rdi = 0.5 * repro * M.erepro[s] / w[wMinIdx[s]];
        r.rdi[s] = rdi;
        r.rdd[s] = rdi / (1 + rdi / M.R_max[s]);
      }
      for (let k = 0; k < NF; k++) {
        let rm = 0;
        for (let i = 0; i < S; i++) rm += thetaR[i] * r.pred_rate[i * NF + k];
        r.resource_mort[k] = rm;
      }
      return r;
    }

    function step(dt) {
      const n = state.n, n_pp = state.n_pp;
      rates(n, n_pp, state.effort);
      // Resource: exact solution with mortality frozen over the step
      for (let k = 0; k < NF; k++) {
        const mur = rr[k] + r.resource_mort[k];
        const nst = rr[k] * cc[k] / mur;
        const v = nst + (n_pp[k] - nst) * Math.exp(-mur * dt);
        if (Number.isFinite(v)) n_pp[k] = v;
      }
      // Fish: implicit upwind (project_n_no_diffusion)
      for (let s = 0; s < S; s++) {
        const o = s * NW, j0 = wMinIdx[s];
        const g = r.e_growth, mu = r.mort;
        let prev = (n[o + j0] + r.rdd[s] * dt / dw[j0]) /
          (1 + g[o + j0] * dt / dw[j0] + mu[o + j0] * dt);
        n[o + j0] = prev;
        for (let j = j0 + 1; j < NW; j++) {
          const a = -g[o + j - 1] * dt / dw[j];
          const b = 1 + g[o + j] * dt / dw[j] + mu[o + j] * dt;
          prev = (n[o + j] - a * prev) / b;
          n[o + j] = prev;
        }
      }
      state.t += dt;
    }

    // ---- Diagnostics -------------------------------------------------------

    function biomass(s, wFrom = 0) {
      let b = 0;
      for (let j = 0; j < NW; j++) {
        if (w[j] >= wFrom) b += state.n[s * NW + j] * w[j] * dw[j];
      }
      return b;
    }

    // Yield rate (biomass caught per year) for species s at the current state
    function yieldRate(s) {
      let y = 0;
      for (let j = 0; j < NW; j++) {
        const sj = s * NW + j;
        y += state.effort[s] * Fpe[sj] * state.n[sj] * w[j] * dw[j];
      }
      return y;
    }

    // Steady cohort per unit egg density under the current growth and
    // mortality (mizer:::get_steady_state_n with upwind flux). Returns the
    // invasion reproduction number R0 and survival to selected sizes.
    function cohort(s, rr_) {
      const o = s * NW, j0 = wMinIdx[s];
      const g = rr_.e_growth, mu = rr_.mort;
      const nc = new Float64Array(NW);
      nc[j0] = 1;
      for (let j = j0 + 1; j < wTop[s]; j++) {
        nc[j] = g[o + j - 1] * nc[j - 1] / (g[o + j] + mu[o + j] * dw[j]);
      }
      const recruits = g[o + j0] + mu[o + j0] * dw[j0];
      let repro = 0;
      for (let j = 0; j < NW; j++) repro += nc[j] * rr_.e_repro[o + j] * dw[j];
      const R0 = 0.5 * M.erepro[s] * repro / M.w_min[s] / recruits;
      const survivalTo = (wt) => {
        let j = 0;
        while (j < NW - 1 && w[j] < wt) j++;
        return g[o + j] * nc[j] / recruits;
      };
      return { R0, survivalTo, n: nc };
    }

    function diagnostics() {
      const rr_ = rates(state.n, state.n_pp, state.effort);
      return { rates: rr_, cohort: [cohort(0, rr_), cohort(1, rr_)] };
    }

    // ---- State manipulation ------------------------------------------------

    function loadState(name) {
      const st = M.states[name];
      for (let s = 0; s < S; s++) {
        for (let j = 0; j < NW; j++) state.n[s * NW + j] = st.n[s][j];
      }
      state.n_pp.set(st.n_pp);
      state.t = 0;
    }

    // Add `factor` times the State RFP spectrum of species s (a restocking seed)
    function seed(s, factor) {
      const ref = M.states.State_RFP.n[s];
      for (let j = 0; j < NW; j++) state.n[s * NW + j] += factor * ref[j];
    }

    function remove(s) {
      for (let j = 0; j < NW; j++) state.n[s * NW + j] = 0;
    }

    function snapshot() {
      return { n: Float64Array.from(state.n), n_pp: Float64Array.from(state.n_pp),
               effort: Float64Array.from(state.effort), t: state.t };
    }

    function restore(snap) {
      state.n.set(snap.n); state.n_pp.set(snap.n_pp);
      state.effort.set(snap.effort); state.t = snap.t;
    }

    return {
      model: M, state, step, rates, biomass, yieldRate, cohort, diagnostics,
      loadState, seed, remove, snapshot, restore,
      dims: { S, NW, NF, off },
    };
  }

  return { createEngine };
});
