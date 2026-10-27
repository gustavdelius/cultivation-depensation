# Project notes: run-simulation

- The model uses physical size diffusion (`D_ext`, c = dσ²/d ln w = 0.1) with
  `second_order_w(params) <- "centred"` and `project(method = "tr_bdf2")`.
  The first-order upwind scheme adds numerical diffusion equivalent to
  c ≈ Δln w, larger than the physical diffusion on an 80-class grid.
- Second-order convergence needs every input to be smooth on the grid: the
  sigmoid (not knife-edge) selectivity, a smooth resource taper instead of a
  sharp `w_pp_cutoff`, and a non-steep maturity ogive. Before mizer#619 was
  fixed (Oct 2026, mizer 3.4.0.9000), predation mortality was only first
  order whatever the inputs.
- Diffusion continues above w_inf, where mean growth is zero. Without the
  senescence term in `build_cultivation_model()` adults spread to several
  times w_inf and results depend on `w_max`. Check boundary independence
  (vary `w_max_factor`) after any change to diffusion or adult mortality.
- The resource grid must extend to `min_w_pp = 1e-9` because encounter and
  diffusion are circular convolutions (mizer#618).
