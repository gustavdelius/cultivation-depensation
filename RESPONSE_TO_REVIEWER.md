# Response to the reproducibility review

We thank the reviewer for an independent reconstruction that was careful enough
to find two real problems: the robustness variants were not calibrated
consistently, and two numbers in the manuscript printed as `NA`. Both are
fixed. The reference model and every result that depends only on it are
unchanged. The changes are in `cultivation_depensation.R` (new
`recalibrate()`, extended `variant_summary()`, new `run_robustness()` and a
`--robustness` mode) and `paper.Rmd`.

The original author files (`cultivation_depensation.R`,
`cultivation_depensation_params.rds`, `results/cultivation_results.rds`) were
the starting point, and their earlier versions remain in git.
`cultivation_depensation_params.rds` was not changed.

## R1 Calibration procedure and the reference model

**Evidence.** The reference model is built by `build_cultivation_model()`:

1. `newMultispeciesParams()` with `biomass_observed = c(100, 600)`;
2. `steadySingleSpecies()`;
3. two rounds of `matchBiomasses()` followed by
   `tuneSteadyState(dt = 0.2, info_level = 0)` (projection solver, default
   tolerances, `R_max = Inf`);
4. `reproduction_level(params) <- c(Forage = 0.65, Predator = 0.78)`;
5. `findSteadyState(solver = "newton")`, which changes no parameter.

There is no stopping test. Two rounds stop short of the targets (93.8/584.3 g),
and the Newton step in 5 moves the state slightly, so the levels become
0.6511/0.7801.

**Answer to R1.2.** The shortfall was accepted calibration error, not
intentional. The published reference model is the output of this recipe,
saved as `cultivation_depensation_params.rds`. A clean rerun reproduces it
exactly: the reference row of the new robustness run has relative biomass
error 1.5e-10 against the saved model and needed 0 extra calibration rounds.

**Sensitivity.** Calibrating the reference model fully to the nominal targets
is now a row of Table 3. It reproduces the reviewer's numbers:

| | Our value | Reviewer |
|---|---:|---:|
| R0_Pred(RF) | 0.220 | 0.2191 |
| F_Forage,inv | 0.340 | 0.3408 |
| F_crit | 0.153 | — |
| F_Forage,rec | 0.394 | — |

**Manuscript changes.**
- Section 2.2 now gives the exact recipe and says that the achieved state,
  not the targets, defines the reference model.
- A new "Calibration" paragraph in Robustness, and the Limitations section,
  report the sensitivity above.
- Reproducibility describes the saved per-variant models and tables.

**R1.3 (per-variant table and models).** Every variant is saved with
`saveParams()`, including resource capacity and numerical settings:
calibrated variants in `results/models/`, recipe-only variants in
`results/models/recipe/`. `results/robustness.csv`,
`results/diffusion_scan.csv` and `results/cannibalism_scan.csv` give for each
variant:
- achieved State RFP biomasses;
- ε, R_max and L_R for each species;
- the number of calibration rounds and the final relative biomass error;
- the biomass drift over a 100-year projection.

**R1.4 (mizer build).** The installed mizer 3.4.0.9001 was built on 28 Oct
2026 from commit `e1e85b0` of sizespectrum/mizer. That commit is on the
unmerged branch of PR #627 (FFT zero-padding, issue #618), not on master,
where the version is still 3.4.0.9000. The Reproducibility section now says
so and gives `remotes::install_github("sizespectrum/mizer@e1e85b0")`.
`results/session_info.txt` and `res$robustness_session` record the session.

**Status: resolved.** The commit was identified from the mizer checkout's
reflog (it was checked out at 11:02 and the package was built at 11:22, with
a clean tree). The installed package itself does not record its commit.

## R2 Predator-cannibalism threshold

**Evidence.**

- **Control.** Changing θ_Pred,Pred alone in the reference State RF, without
  recalibrating, leaves R0_Pred(RF) = 0.27552 at θ = 0.15, 0.36 and 0.95. The
  predator is absent from RF, so cannibalism does not enter the invasion
  calculation. The whole effect in the scan therefore comes from
  recalibration.
- **Original protocol.** The fixed recipe alone falls increasingly short of
  the targets as cannibalism increases. State RFP holds 271–584 g of
  predators and 54–94 g of forage fish across the scan, which also changes
  the resource capacity and State RF. Rerunning that protocol reproduces the
  published thresholds exactly: 0.149 (θ_PP) and 0.118 (θ_FF).
- **Converged protocol.** Each variant is calibrated to the reference State
  RFP: biomass error below 1e-6, then the reference reproduction levels are
  imposed and verified to 1e-5. The thresholds are found by `uniroot` on
  separately calibrated models (tolerance 1e-4).

| | Converged (new) | Reviewer | Recipe only (old) |
|---|---:|---:|---:|
| θ_For,For | 0.1203 | 0.12030 | 0.118 |
| θ_Pred,Pred | 0.3596 | 0.35964 | 0.149 |

- **Mechanism.** In the converged scan State RF barely changes
  (2169–2435 g of forage fish). What changes is the predator's reproductive
  efficiency, which the calibration raises from 3.0e-5 at θ_PP = 0.05 to
  1.5e-4 at θ_PP = 0.4 to offset cannibalism in RFP. An invader in RF keeps
  that efficiency without suffering cannibalism.

**Manuscript changes.**
- The abstract and results now give θ_Pred,Pred = 0.360. The forage value
  goes from 0.118 to 0.120.
- The predator-cannibalism paragraph is rewritten around the control and the
  mechanism, and reports both protocols.
- Figure 7 shows both protocols: the converged scan in orange, the recipe
  only in grey.
- The old parenthesis about the calibration missing its target is removed.

**Status: resolved.** The 0.149 value is explained, reproduced and replaced.
No parameters were tuned to match any value.

## R3 Diffusion sensitivity

**Evidence.** The published diffusion variants had the same problem as R2:
State RFP biomasses differed from the reference (for example, 589.6 g of
predators at c = 0.02). With converged calibration, the same transport
schemes (van Leer at c = 0.02, centred otherwise) and the reference
reproduction levels, the new values match the reviewer's:

| c | R0 (ours / reviewer) | F_crit | F_Forage,inv | F_Forage,rec |
|---:|---|---|---|---|
| 0.02 | 0.2141 / 0.2139 | 0.1584 / 0.1585 | 0.2920 / 0.2922 | 0.3270 / 0.3268 |
| 0.05 | 0.2393 / 0.2391 | 0.1628 / 0.1632 | 0.2921 / 0.2922 | 0.3332 / 0.3332 |
| 0.10 | 0.2755 / 0.2753 | 0.1697 / 0.1696 | 0.2880 / 0.2881 | 0.3370 / 0.3373 |
| 0.20 | 0.3444 / 0.3442 | 0.1834 / 0.1837 | 0.2712 / 0.2713 | 0.3342 / 0.3346 |
| 0.30 | 0.4131 / 0.4129 | 0.1991 / 0.1995 | 0.2489 / 0.2491 | 0.3233 / 0.3236 |

The remaining differences are within the bisection tolerance of the
thresholds (0.001).

**Manuscript changes.** Figure 6 and the quoted ranges are regenerated from
these values, and the point-wise values are in `results/diffusion_scan.csv`.
The ranges in the text are now:
- R0: 0.21–0.41;
- F_crit: 0.158–0.199;
- F_Forage,inv: 0.249–0.292;
- F_Forage,rec: 0.323–0.337.

The caption now names the van Leer limiter at c = 0.02.

**Status: resolved.**

## R4 Convergence diagnostic at the upper size boundary

**Evidence.** We confirm the behaviour in our own reference model. The
predator's last size class (320) lies outside the simulated support
(`support_top_idx` = 319) and holds zero density, yet `getSteadyResidual()`
assigns it 1.6e-8/yr. That is negligible here (`isSteady()` is TRUE), but it
can dominate at other parameter values, as in the reviewer's c = 0.2 example.

The analysis never relied on `getSteadyResidual()`. Steady states are checked
by the drift of total biomass over a 100-year TR-BDF2 projection (below
2.2e-9/yr for every state and every calibrated variant; see the `drift_RFP`
column of the CSV files).

**Manuscript change.** Section 2.3 now states this diagnostic and why
`getSteadyResidual()` is not used.

**Status: resolved for the paper.** The behaviour itself is a mizer issue.
The reviewer's draft is a reasonable basis for reporting it at
sizespectrum/mizer, but it has not been posted.

## R5 Missing values and unverified items

- **The two `NA` values.** These came from a bug in `paper.Rmd`, not in the
  calculation: it indexed `res$R0_decomposition` by `"growthRFP_mortRF"`,
  but the stored names carry a `.Predator` suffix. The manuscript now prints
  0.136 (RFP growth, RF mortality) and 11.8 (RF growth, RFP mortality),
  matching the reviewer's 0.135526 and 11.780840.
- **The substitutions.** `R0_mix()` takes g_Pred and E_r ψ from the
  "growth" state, μ_Pred (total mortality) from the other, and diffusion from
  the growth state (it does not depend on the state). The text now says
  exactly this.
- **Check on the decomposition.** Without substitution the same calculation
  returns 4.547 in RFP, equal to 1/(1 − L_R), and 0.276 in RF. Both are now
  quoted.
- **Survival definition.** Survival to maturation is now defined in the text
  as exp(−∫ μ/g dw) from w_min to w_mat along the mean growth trajectory,
  summed over size classes whose lower edge is ≤ w_mat and ignoring size
  diffusion.
- **Grid-resolution, first-order and resource-taper variants.** These are
  recomputed with the converged calibration, and their models and values are
  saved. Their numbers change slightly: the first-order R0 goes from 0.476 to
  0.493 and the w_R = 1 g F_crit from 0.652 to 0.467. The conclusions drawn
  from them in the text still hold.

**Status: resolved.** We have not claimed that the paper as a whole has been
independently reproduced. Verifying these variants independently remains
open to the reviewer.
