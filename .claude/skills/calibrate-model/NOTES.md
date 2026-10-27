# Project notes: calibrate-model

- With size diffusion, `getSteadyResidual()` and the `residual` in
  `attr(params, "convergence")` report per-size-class relative rates that the
  nearly empty classes near w_inf and w_max dominate (0.05 for a converged
  state). Judge convergence by the biomass drift (`biomass_drift()` in
  `cultivation_depensation.R`) instead.
- `mizer:::get_steady_state_n()` stops silently at 500 iterations; the
  scheme-consistent R0 in the script passes `max_iterations = 20000`.
