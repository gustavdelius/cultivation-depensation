# Project notes: calibrate-model

- With size diffusion, `getSteadyResidual()` and the `residual` in
  `attr(params, "convergence")` report per-size-class relative rates that the
  nearly empty classes near w_inf and w_max dominate (0.05 for a converged
  state). Judge convergence by the biomass drift (`biomass_drift()` in
  `cultivation_depensation.R`) instead.
- `mizer:::get_steady_state_n()` stops silently at 500 iterations; the
  scheme-consistent R0 in the script passes `max_iterations = 20000`.
- A fixed number of `matchBiomasses()`/`tuneSteadyState()` rounds is not a
  calibration: the reference recipe (two rounds) stops 6%/3% short of its
  100/600 g targets, and much further for variants with strong cannibalism,
  which made variant comparisons confound the varied parameter with a
  different State RFP (predator-cannibalism threshold 0.149 vs 0.360).
  `recalibrate()` in `cultivation_depensation.R` alternates `matchBiomasses()`
  with `tuneSteadyState(solver = "newton")` to a relative error of 1e-6
  (4-45 rounds), then sets `reproduction_level` and verifies both. Biomass
  error shrinks by a constant factor (~0.4-0.6) per round.
- `tuneSteadyState()` (default `preserve = "reproduction_level"`) still lets
  the reproduction levels drift during these rounds; impose them only after
  the biomass calibration and check them afterwards.
