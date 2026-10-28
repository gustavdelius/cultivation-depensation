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
- `min_w_pp = 1e-9` puts the larval feeding kernels on the grid to more than
  six standard deviations. It was first chosen to keep clear of mizer#618
  (circular FFT convolutions wrapping around the size grid), which mizer
  3.4.0.9001 fixed. Re-running with the fix (Oct 2026) changed the reference
  model only at round-off; the "Original model" variant (`min_w_pp = 1e-6`)
  moved by 2-4%. Do not shorten the grid to `1e-6`: that cuts off the lower
  tail of the larval kernels.
- The full run takes about 110 min with `MC_CORES=8` on the user's laptop
  (4 physical cores). Calibration, the saddle, the invasion and management
  runs, the F_Forage continuation and the F_forage_rec bisection are serial.
  `mclapply` preschedules jobs, so the regime map ends with one worker still
  running the slow jobs near the cusp.
- After rerunning, render the paper and look at the figures, not only the
  results file: before Oct 2026 the committed figures came from older runs
  than the committed results (old state names, a different invasion decay).
