# Project notes: change-parameters

- To vary the strength of competition for the resource without recalibrating,
  scale the regeneration rate with `resource_rate(p, balance = TRUE) <- ...`
  (`with_resource_rate()` in `cultivation_depensation.R`). The capacity is
  rebalanced so State RFP and every fish parameter stay exactly as they were;
  at r0 = 10 this reproduces the reference thresholds to all digits.
  `resource_dynamics(p) <- "resource_constant"` cannot stand in for the
  no-competition limit here: `findSteadyState(solver = "newton")`, used by
  `polish()`, refuses non-semichemostat dynamics. Use a large r0 (1000).
- `--figures` redraws `plots/predator_exploitation_bifurcation.png` with
  different bytes from the committed file even when nothing changed; restore
  it with `git checkout` unless the branches really changed.
