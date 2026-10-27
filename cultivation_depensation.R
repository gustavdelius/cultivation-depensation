# Cultivation-Depensation in a 3-Taxa Mizer Model
# Community:
#   1. Plankton resource spectrum (semichemostat dynamics)
#   2. Forage fish spectrum (small maximum size, w_inf = 20 g)
#   3. Predator size spectrum (large maximum size, w_inf = 2000 g)
#
# Run with `Rscript cultivation_depensation.R`. It writes the figures to plots/,
# every number quoted in paper.Rmd to results/cultivation_results.rds and the
# calibrated model to cultivation_depensation_params.rds. The full run takes
# roughly half an hour on four cores.

suppressPackageStartupMessages({
  library(mizer)
  library(ggplot2)
  library(patchwork)
})
if (packageVersion("mizer") < "3.4.0") {
  stop("This script needs mizer >= 3.4.0; the saved model cannot be read by older versions.")
}

#' Lifetime reproductive output (R0) of a rare individual in a resident environment
#'
#' Holding growth, mortality and reproductive investment at their values in
#' `params`, solves the steady-state transport equation of the model's own
#' size discretisation for the density produced by one recruit, and returns the
#' density-independent egg production of that density. The recruitment flux
#' needed to sustain the density is read from the balance of the first size
#' class, so the result is consistent with whichever transport scheme the model
#' uses. A rare species can invade a resident steady state iff R0 > 1, and for
#' a resident species at a Beverton-Holt steady state R0 = 1 / (1 - reproduction_level).
get_R0 <- function(params, effort = initial_effort(params),
                   species = species_params(params)$species) {
  p <- params
  initial_effort(p) <- effort
  R0_from_rates(p, getEGrowth(p), getMort(p), getERepro(p), getDiffusion(p), species)
}

#' R0 for given growth, mortality, reproductive investment and diffusion rates
R0_from_rates <- function(p, g, mort, e_repro, D, species) {
  sp_names <- species_params(p)$species
  flux_limiter <- mizer:::flux_limiter_scheme(p)
  R0 <- setNames(numeric(length(species)), species)
  for (sp in species) {
    N0 <- setNames(numeric(length(sp_names)), sp_names)
    N0[sp] <- 1
    # The flux limiter makes the solve iterative; the default iteration cap is
    # too small for it to converge here.
    n <- mizer:::get_steady_state_n(p, g, mort, D, N0, max_iterations = 20000)
    coefs <- mizer:::get_transport_coefs(p, n, g, mort, dt = 1,
                                         recruitment_flux = rep(0, length(sp_names)),
                                         d = D, flux_limiter = flux_limiter)
    j <- p@w_min_idx[sp]
    R_req <- dw(p)[j] * ((coefs$b[sp, j] - 1) * n[sp, j] + coefs$c[sp, j] * n[sp, j + 1])
    rdi <- mizerRDI(p, n = n, n_pp = initialNResource(p), n_other = initialNOther(p),
                    e_growth = g, mort = mort, e_repro = e_repro, diffusion = D)
    R0[sp] <- unname(rdi[sp] / R_req)
  }
  R0
}

#' Build and calibrate the cultivation-depensation MizerParams at Coexistence State A
#'
#' The defaults give the reference model. The arguments exist for the
#' robustness checks: grid resolution, intraspecific predation, the upper size
#' of the resource spectrum and the strength of size diffusion.
#'
#' `c_diff` sets the physical size diffusion D(w) = D_ext w^(n + 1), which
#' represents individual variability in growth. mizer's transport equation
#' contains the diffusion as (1/2) d^2(D N) / dw^2, so the variance of body
#' mass in a cohort grows at rate D and the variance of ln w grows by
#' d sigma^2 / d ln w = D / (g w) per unit of ln w grown. With juvenile growth
#' g(w) ~ (alpha f0 h - ks) w^n, D_ext is chosen to make this equal to `c_diff`. The size axis is discretised with mizer's second-order scheme so
#' that numerical diffusion stays small compared with this physical diffusion.
#' The resource grid starts far below the larval prey sizes so that the FFT
#' convolutions do not wrap around (sizespectrum/mizer#618).
#'
#' Fishing selectivity is logistic in length, with 50% selection at 1 g
#' (forage) and 400 g (predator) and 25% selection at 90% of those lengths.
#' Likewise 25% of fish are mature at `ogive_l25` times the length at 50%
#' maturity; the default 0.9 gives a 25-75% maturation range of about 20% of
#' that length, whereas mizer's default ogive (`ogive_l25 = NA`) is much steeper.
#'
#' The resource capacity tapers off smoothly around `w_pp_cutoff`,
#' kappa w^-lambda / (1 + (w / w_pp_cutoff)^taper); a sharp cutoff
#' (`taper = NA`) puts a discontinuity into the prey spectrum whose position
#' relative to the size grid makes the feeding integrals only first-order
#' accurate.
#'
#' Size diffusion continues where mean growth stops, so without further
#' mortality adults spread far beyond w_inf and the results depend on where
#' the grid ends (`w_max = w_max_factor * w_inf`). A senescence mortality
#' z_sen (w / w_inf)^k_sen removes that tail.
#'
#' `flux = "upwind"` switches to mizer's original first-order scheme; together
#' with `no_w = 80`, `c_diff = 0`, `min_w_pp = 1e-6`,
#' `selectivity = "knife_edge"`, `ogive_l25 = NA`, `taper = NA`,
#' `w_max_factor = 1` and `z_sen = 0` it reproduces the model of the first
#' version of this study.
build_cultivation_model <- function(no_w = 320, theta_FF = 0.05, theta_PP = 0.05,
                                    w_pp_cutoff = 10, c_diff = 0.1,
                                    flux = "centred", min_w_pp = 1e-9,
                                    selectivity = "sigmoid", ogive_l25 = 0.9,
                                    taper = 4, w_max_factor = 2,
                                    z_sen = 0.5, k_sen = 10) {
  sp <- data.frame(
    species          = c("Forage", "Predator"),
    w_inf            = c(20, 2000),
    w_max            = w_max_factor * c(20, 2000),
    w_mat            = c(4, 250),
    w_min            = c(0.001, 0.001),
    beta             = c(200, 100),
    sigma            = c(1.3, 1.5),
    k_vb             = c(1.0, 0.5),
    f0               = c(0.6, 0.6),
    z0               = c(0.30, 0.20),
    z_ext            = c(0.50, 1.20),
    a                = 0.01,  # length-weight relation w = a l^b (w in g, l in cm)
    b                = 3,
    biomass_observed = c(100, 600)
  )
  if (!is.na(ogive_l25)) sp$w_mat25 <- ogive_l25^sp$b * sp$w_mat

  # Interaction matrix:
  # - Predator preys strongly on Forage (1.00): top-down cultivation
  # - Forage preys on Predator larvae/early juveniles (0.95): depensatory bottleneck
  # - Weak self-interaction (0.05) avoids single-species cannibalism oscillations
  inter <- matrix(
    c(theta_FF, 0.95,
      1.00,     theta_PP),
    nrow = 2, byrow = TRUE,
    dimnames = list(sp$species, sp$species)
  )

  # Separate species-specific gears with catchability = 1, so that effort equals
  # the fishing mortality on fully selected fish (1/year)
  w50 <- c(1.0, 400.0)
  gp <- data.frame(
    gear            = c("Forage", "Predator"),
    species         = c("Forage", "Predator"),
    catchability    = c(1.0, 1.0)
  )
  if (selectivity == "knife_edge") {
    gp$sel_func <- "knife_edge"
    gp$knife_edge_size <- w50
  } else {
    gp$sel_func <- "sigmoid_length"
    gp$l50 <- (w50 / sp$a)^(1 / sp$b)
    gp$l25 <- 0.9 * gp$l50
  }

  params <- newMultispeciesParams(
    sp,
    interaction      = inter,
    gear_params      = gp,
    initial_effort   = c(Forage = 0, Predator = 0),
    no_w             = no_w,
    min_w_pp         = min_w_pp,
    kappa            = 100,
    w_pp_cutoff      = if (is.na(taper)) w_pp_cutoff else 100 * w_pp_cutoff,
    resource_rate    = 10,
    second_order_w   = flux != "upwind",
    info_level       = 0
  )
  if (flux != "upwind") second_order_w(params) <- flux
  if (!is.na(taper)) {
    # Capacity values are bin averages; apply the taper at the bin centres
    wf <- w_full(params)
    w_centre <- wf * exp(log(wf[2] / wf[1]) / 2)
    cap <- resource_capacity(params) / (1 + (w_centre / w_pp_cutoff)^taper)
    params <- setResource(params, resource_capacity = cap, balance = FALSE)
    initialNResource(params) <- cap
  }
  if (z_sen > 0) {
    # Senescence: without it the size diffusion, which continues where mean
    # growth stops, spreads adults far beyond w_inf and the results depend on
    # w_max. Bin averages of (w / w_inf)^k_sen on a second-order grid.
    w_l <- w(params); w_r <- w_l + dw(params)
    sen <- outer(z_sen / species_params(params)$w_inf^k_sen,
                 if (flux == "upwind") w_l^k_sen else
                   (w_r^(k_sen + 1) - w_l^(k_sen + 1)) / ((k_sen + 1) * dw(params)))
    ext_mort(params) <- ext_mort(params) + sen
  }
  sp_p <- species_params(params)
  given_species_params(params)$D_ext <-
    c_diff * (sp_p$alpha * sp_p$f0 * sp_p$h - sp_p$ks)

  params <- setColours(params, c(Forage = "#0072B2", Predator = "#D55E00"))

  # Calibrate to Coexistence State A
  params <- steadySingleSpecies(params)
  params <- matchBiomasses(params)
  params <- tuneSteadyState(params, dt = 0.2, info_level = 0)
  params <- matchBiomasses(params)
  params <- tuneSteadyState(params, dt = 0.2, info_level = 0)

  # Set Beverton-Holt recruitment density dependence at State A
  reproduction_level(params) <- c(Forage = 0.65, Predator = 0.78)

  # Solve for fixed point A without altering parameters
  params_A <- findSteadyState(params, solver = "newton", info_level = 0)
  params_A <- setMetadata(
    params_A,
    title = "Cultivation-Depensation 3-Taxa Mizer Model",
    description = paste(
      "A minimal 3-taxa size-spectrum model (Plankton, Forage fish, Predator)",
      "exhibiting stable coexistence (State A), a plankton-predator state (State B)",
      "invasible by forage fish, and a plankton-forage-fish state (State C)",
      "uninvasible by the predator without forage-fish harvesting."
    )
  )
  params_A
}

#' Project towards a fixed point and finish with the Newton solver
polish <- function(params, effort = initial_effort(params)) {
  p <- finalParams(project(params, effort = effort, t_max = 400, dt = 0.5,
                           progress_bar = FALSE))
  findSteadyState(p, solver = "newton", effort = effort, info_level = 0)
}

#' Derive the four unfished fixed points (State 0, State A, State B, State C)
#' sharing identical model parameters
get_fixed_points <- function(params_A) {
  # State B: Plankton + Predator subsystem (Forage absent)
  params_B <- params_A
  initialN(params_B)["Forage", ] <- 0
  params_B <- polish(params_B)

  # State C: Plankton + Forage subsystem (Predator absent)
  params_C <- params_A
  initialN(params_C)["Predator", ] <- 0
  params_C <- polish(params_C)

  # State 0: Plankton only (both fish absent, resource at carrying capacity)
  params_0 <- params_A
  initialN(params_0)[] <- 0
  initialNResource(params_0) <- resource_capacity(params_0)

  list(State_0 = params_0, State_A = params_A, State_B = params_B, State_C = params_C)
}

# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------

# Predator biomass [g] separating the abundant (cultivated) states from the
# depleted ones (State C or a low-abundance coexistence state).
B_HI <- 100

#' Project to a steady state
#'
#' Steady states of mizer's scheme do not depend on `dt`, so the long settling
#' runs use `dt = 0.5`. The run is extended until every species' biomass changes
#' by less than `tol` per year, or is below 1e-4 g and still declining. The
#' tolerance is strict because close to a fold the dynamics pass slowly through
#' the remnant of the vanished equilibrium.
settle <- function(params, effort, t_max = 500, max_rounds = 10, tol = 1e-5) {
  p <- params
  for (k in seq_len(max_rounds)) {
    sim <- project(p, effort = effort, t_max = t_max, dt = 0.5,
                   t_save = t_max / 10, progress_bar = FALSE)
    p <- finalParams(sim)
    b <- getBiomass(sim)
    nt <- nrow(b)
    rate <- log(b[nt, ] / b[nt - 1, ]) / (t_max / 10)
    if (all(is.na(rate) | abs(rate) < tol | (b[nt, ] < 1e-4 & rate < 0))) break
  }
  p
}

pred_biomass <- function(params) unname(getBiomass(params)["Predator"])

#' Largest relative rate of biomass change [1/year] over a 100-year projection
#'
#' Used to confirm that a state is steady. (`getSteadyResidual()` reports
#' relative rates in individual size classes, which are dominated by nearly
#' empty classes at the top of the size range once there is size diffusion.)
biomass_drift <- function(params, effort) {
  sim <- project(params, effort = effort, t_max = 100, dt = 0.1, t_save = 100,
                 method = "tr_bdf2", progress_bar = FALSE)
  b <- getBiomass(sim)
  present <- b[1, ] > 0
  max(abs(log(b[2, present] / b[1, present])) / 100)
}

#' Eigenvalues of the linearised dynamics, without unreachable size classes
#'
#' Species named in `invaders` are added at trace density, with the size
#' distribution they have in `ref`, so that their invasion exponent is part of
#' the spectrum. Size classes above the largest size a species can grow to hold
#' no individuals; each contributes an eigenvalue equal to minus its mortality
#' rate, which describes no ecological dynamics, so these are removed.
dynamic_eigenvalues <- function(params, effort, invaders = character(),
                                ref = params) {
  p <- params
  for (sp in invaders) initialN(p)[sp, ] <- initialN(ref)[sp, ] * 1e-12
  initial_effort(p) <- effort
  ev <- suppressWarnings(getStability(p, effort = effort))$eigenvalues
  n <- initialN(p)
  mort <- getMort(p, effort = effort)
  empty <- (n == 0) & (rowSums(n) > 0)
  keep <- rep(TRUE, length(ev))
  for (v in -mort[empty]) {
    d <- abs(ev - v)
    d[!keep] <- Inf
    j <- which.min(d)
    if (d[j] < 1e-6 * max(1, abs(v))) keep[j] <- FALSE
  }
  ev[keep]
}

leading_eigenvalue <- function(...) max(Re(dynamic_eigenvalues(...)))

#' Unstable interior equilibrium separating States A and C
#'
#' Edge tracking: bisect on the size of a predator seed added to State C to find
#' the basin boundary, follow a trajectory started on it to where it lingers
#' near the saddle, and polish that state with the Newton solver.
find_saddle <- function(params_A, params_C, effort = c(Forage = 0, Predator = 0)) {
  outcome <- function(s) {
    p <- params_C
    initialN(p)["Predator", ] <- s * initialN(params_A)["Predator", ]
    sim <- project(p, effort = effort, t_max = 400, dt = 0.5, t_save = 400,
                   progress_bar = FALSE)
    getBiomass(sim)[2, "Predator"] > B_HI
  }
  lo <- 0; hi <- 1
  for (k in 1:30) {
    m <- (lo + hi) / 2
    if (outcome(m)) hi <- m else lo <- m
  }
  p <- params_C
  initialN(p)["Predator", ] <- lo * initialN(params_A)["Predator", ]
  sim <- project(p, effort = effort, t_max = 400, dt = 0.5, t_save = 1,
                 progress_bar = FALSE)
  b <- getBiomass(sim)[, "Predator"]
  drift <- abs(diff(log(b)))
  k <- which.min(drift[-(1:5)]) + 6
  initialN(p) <- N(sim)[k, , ]
  initialNResource(p) <- NResource(sim)[k, ]
  saddle <- findSteadyState(p, solver = "newton", effort = effort, info_level = 0)
  ev <- dynamic_eigenvalues(saddle, effort)
  list(params = saddle, seed_scale = lo, biomass = getBiomass(saddle),
       convergence = attr(saddle, "convergence"),
       leading = max(Re(ev)), n_unstable = sum(Re(ev) > 0))
}

#' R0 of a rare species using growth from one environment and mortality from another
R0_mix <- function(p_growth, p_mort, sp = "Predator") {
  R0_from_rates(p_growth, getEGrowth(p_growth), getMort(p_mort), getERepro(p_growth),
                getDiffusion(p_growth), sp)
}

#' Follow the stable states along a sequence of predator fishing mortalities
#'
#' Each step starts from the state reached at the previous one. With `seed`
#' the predator is topped up to at least that abundance before each step, so
#' that it can invade wherever its invasion exponent is positive.
trace_branch <- function(start, Ff, Fp_seq, seed = NULL, t_max = 500) {
  p <- start
  out <- vector("list", length(Fp_seq))
  for (i in seq_along(Fp_seq)) {
    if (!is.null(seed)) initialN(p)["Predator", ] <- pmax(initialN(p)["Predator", ], seed)
    p <- settle(p, c(Forage = Ff, Predator = Fp_seq[i]), t_max = t_max)
    b <- getBiomass(p)
    out[[i]] <- data.frame(F_Forage = Ff, F_Pred = Fp_seq[i],
                           Forage = unname(b["Forage"]), Predator = unname(b["Predator"]))
  }
  do.call(rbind, out)
}

#' Locate the regime boundaries along the F_Pred axis at fixed F_Forage
#'
#' * `F_inv`: predator invasion boundary of State C, R0_Predator(State C) = 1.
#' * `F_up`: largest F_Pred at which the abundant coexistence state persists
#'   (continuation from State A, refined by bisection). `lead_up` is the
#'   leading eigenvalue just below it; a value near zero identifies a fold.
#' * `F_low`: smallest F_Pred at which the depleted state (State C, or a
#'   low-abundance coexistence state `A_low` below `F_inv`) persists
#'   (continuation from State C with a predator seed, refined by bisection).
map_regime <- function(params_A, params_C, Ff, step = 0.01, tol = 0.001) {
  eff <- function(Fp) c(Forage = Ff, Predator = Fp)
  seed <- initialN(params_A)["Predator", ] * 1e-6

  # Invasion boundary of State C (State C does not depend on F_Pred)
  pC <- settle(params_C, eff(0))
  R0_C0 <- unname(get_R0(pC, effort = eff(0), species = "Predator")["Predator"])
  F_inv <- if (R0_C0 > 1) {
    uniroot(function(fp) get_R0(pC, effort = eff(fp), species = "Predator")["Predator"] - 1, c(0, 2))$root
  } else NA_real_

  # Upper limit of the abundant state
  p <- settle(params_A, eff(0)); Fp <- 0
  while (Fp < 2) {
    q <- settle(p, eff(Fp + step))
    if (pred_biomass(q) < B_HI) break
    p <- q; Fp <- Fp + step
  }
  lo <- Fp; hi <- Fp + step
  while (hi - lo > tol) {
    m <- (lo + hi) / 2
    q <- settle(p, eff(m))
    if (pred_biomass(q) > B_HI) { lo <- m; p <- q } else hi <- m
  }
  F_up <- (lo + hi) / 2
  lead_up <- leading_eigenvalue(p, eff(lo))
  B_up <- pred_biomass(p)

  # Lower limit of the depleted state
  A_low <- FALSE
  if (is.na(F_inv)) {
    F_low <- 0
  } else {
    p <- pC; last_low <- NULL; F_hi_found <- NA
    for (Fp in seq(F_inv + step, 0, by = -step)) {
      initialN(p)["Predator", ] <- pmax(initialN(p)["Predator", ], seed)
      q <- settle(p, eff(Fp), t_max = 1000)
      if (pred_biomass(q) > B_HI) { F_hi_found <- Fp; break }
      last_low <- list(p = q, Fp = Fp); p <- q
    }
    if (is.na(F_hi_found)) {
      F_low <- 0
    } else if (is.null(last_low)) {
      F_low <- F_hi_found
    } else {
      lo <- F_hi_found; hi <- last_low$Fp; p <- last_low$p
      while (hi - lo > tol) {
        m <- (lo + hi) / 2
        initialN(p)["Predator", ] <- pmax(initialN(p)["Predator", ], seed)
        q <- settle(p, eff(m), t_max = 1000)
        if (pred_biomass(q) > B_HI) lo <- m else { hi <- m; p <- q }
      }
      F_low <- (lo + hi) / 2
    }
    # The depleted state persisting clearly below F_inv is a low-abundance
    # coexistence state: the invasion boundary is then supercritical. A margin
    # of 0.01 allows for the slow growth of the seed close to F_inv.
    A_low <- F_low < F_inv - 0.01
    if (!A_low) F_low <- F_inv
  }
  data.frame(F_Forage = Ff, R0_C = R0_C0, F_inv = F_inv, F_up = F_up,
             lead_up = lead_up, B_up = B_up, F_low = F_low, A_low = A_low)
}

#' Regime at each point of a (F_Forage, F_Pred) grid from the mapped boundaries
classify_regimes <- function(bounds, Fp_grid) {
  out <- lapply(seq_len(nrow(bounds)), function(i) {
    b <- bounds[i, ]
    F_inv <- ifelse(is.na(b$F_inv), -Inf, b$F_inv)
    hysteresis <- b$F_low < b$F_up - 0.003
    regime <- vapply(Fp_grid, function(Fp) {
      high <- Fp < b$F_up
      low <- Fp > b$F_low
      if (high && low && hysteresis) {
        if (Fp > F_inv) "Bistable: A and C" else "Bistable: A and low-predator state"
      } else if (high || (!hysteresis && Fp < F_inv)) {
        "Coexistence only"
      } else if (Fp < F_inv) {
        "Low-predator coexistence only"
      } else {
        "State C only"
      }
    }, character(1))
    data.frame(F_Forage = b$F_Forage, F_Pred = Fp_grid, Regime = regime)
  })
  do.call(rbind, out)
}

#' Collapse threshold F_crit of the abundant state at fixed F_Forage
find_F_crit <- function(params_A, Ff = 0, start = 0.1, step = 0.01, tol = 0.001) {
  eff <- function(Fp) c(Forage = Ff, Predator = Fp)
  p <- settle(params_A, eff(start)); Fp <- start
  while (Fp < 2) {
    q <- settle(p, eff(Fp + step))
    if (pred_biomass(q) < B_HI) break
    p <- q; Fp <- Fp + step
  }
  lo <- Fp; hi <- Fp + step
  while (hi - lo > tol) {
    m <- (lo + hi) / 2
    q <- settle(p, eff(m))
    if (pred_biomass(q) > B_HI) { lo <- m; p <- q } else hi <- m
  }
  (lo + hi) / 2
}

#' Forage-fish mortality at which the predator can invade State C (F_Pred = 0)
find_F_forage_inv <- function(params_C, upper = 1) {
  r0 <- function(Ff) {
    eff <- c(Forage = Ff, Predator = 0)
    get_R0(settle(params_C, eff), effort = eff, species = "Predator") - 1
  }
  if (r0(0) > 0) return(0)
  uniroot(r0, c(0, upper), tol = 1e-4)$root
}

#' Smallest forage-fish mortality at which a predator reintroduced into
#' State C (at 0.5% of its State A abundance, F_Pred = 0) recovers to the
#' abundant state
find_F_forage_rec <- function(params_A, params_C, lower, upper = 1, tol = 0.001) {
  recovers <- function(Ff) {
    eff <- c(Forage = Ff, Predator = 0)
    p <- settle(params_C, eff)
    initialN(p)["Predator", ] <- initialN(params_A)["Predator", ] * 5e-3
    pred_biomass(settle(p, eff, t_max = 1000)) > B_HI
  }
  lo <- lower; hi <- upper
  while (hi - lo > tol) {
    m <- (lo + hi) / 2
    if (recovers(m)) hi <- m else lo <- m
  }
  (lo + hi) / 2
}

#' Summary of a variant of the model, optionally with its fishing thresholds
variant_summary <- function(label, thresholds = FALSE, ...) {
  pA <- suppressMessages(build_cultivation_model(...))
  st <- get_fixed_points(pA)
  diet <- getDiet(pA, proportion = FALSE)
  intake <- colSums(sweep(diet["Predator", , ], 1, initialN(pA)["Predator", ] * dw(pA), "*"))
  out <- data.frame(
    Variant = label,
    B_For_A = unname(getBiomass(st$State_A)["Forage"]),
    B_Pred_A = unname(getBiomass(st$State_A)["Predator"]),
    B_For_C = unname(getBiomass(st$State_C)["Forage"]),
    R0_For_B = unname(get_R0(st$State_B, species = "Forage")),
    R0_Pred_C = unname(get_R0(st$State_C, species = "Predator")),
    Forage_in_diet = unname(intake["Forage"] / sum(intake))
  )
  if (thresholds) {
    out$F_crit <- find_F_crit(pA)
    out$F_forage_inv <- find_F_forage_inv(st$State_C)
    out$F_forage_rec <- find_F_forage_rec(pA, st$State_C, lower = out$F_forage_inv)
  }
  out
}

if (sys.nframe() == 0) {
  dir.create("plots", showWarnings = FALSE)
  dir.create("results", showWarnings = FALSE)
  n_cores <- max(1, min(4, parallel::detectCores()))
  res <- list(mizer_version = as.character(packageVersion("mizer")))

  cat("Building and calibrating cultivation-depensation model...\n")
  params_A <- build_cultivation_model()
  saveParams(params_A, "cultivation_depensation_params.rds")

  states <- get_fixed_points(params_A)
  params_0 <- states$State_0
  params_B <- states$State_B
  params_C <- states$State_C
  F0 <- c(Forage = 0, Predator = 0)

  # -------------------------------------------------------------------------
  # Parameters actually used by the model
  # -------------------------------------------------------------------------
  spp <- species_params(params_A)
  res$species_params <- data.frame(
    species = spp$species, w_inf = spp$w_inf, w_max = spp$w_max,
    w_repro_max = spp$w_repro_max, w_mat = spp$w_mat, w_mat25 = spp$w_mat25, w_min = spp$w_min,
    beta = spp$beta, sigma = spp$sigma, k_vb = spp$k_vb, f0 = spp$f0,
    gamma = spp$gamma, h = spp$h, ks = spp$ks, n = spp$n, p = spp$p, q = spp$q,
    alpha = spp$alpha, z0 = spp$z0, z_ext = spp$z_ext, erepro = spp$erepro,
    R_max = spp$R_max, reproduction_level = unname(reproduction_level(params_A)),
    D_ext = spp$D_ext, l50 = gear_params(params_A)$l50, l25 = gear_params(params_A)$l25
  )
  res$c_diff <- formals(build_cultivation_model)$c_diff
  wf <- w_full(params_A)
  cc <- resource_capacity(params_A)
  rp <- resource_params(params_A)
  bf <- formals(build_cultivation_model)
  # The capacity exceeds the State A resource by the factor 1 + mu_R / r_R;
  # report its maximum over the sizes that fish feed on
  nr <- initialNResource(params_A)
  fed <- nr > 0 & wf >= min(w(params_A)) / 1e4
  res$resource <- list(
    kappa = rp$kappa, lambda = rp$lambda, n = rp$n, r_pp = rp$r_pp,
    w_taper = bf$w_pp_cutoff, taper = bf$taper, max_w_pp = rp$w_pp_cutoff,
    min_w_pp = min(wf),
    capacity_ratio_max = max(cc[fed] / nr[fed]),
    w_at_capacity_ratio_max = wf[fed][which.max(cc[fed] / nr[fed])]
  )
  res$grid <- list(no_w = length(w(params_A)), dx = log(w(params_A)[2] / w(params_A)[1]))
  em <- ext_mort(params_A)
  res$ext_mort_wmin <- em[, params_A@w_min_idx[1]]
  res$senescence <- list(z_sen = bf$z_sen, k_sen = bf$k_sen,
                         w_max_factor = bf$w_max_factor)
  # Share of each species' biomass above w_inf at State A
  bw <- sweep(initialN(params_A), 2, w(params_A) * dw(params_A), "*")
  res$biomass_above_w_inf <- sapply(rownames(bw), function(s)
    sum(bw[s, w(params_A) > species_params(params_A)[s, "w_inf"]]) / sum(bw[s, ]))

  # -------------------------------------------------------------------------
  # Fixed points, R0 and eigenvalues
  # -------------------------------------------------------------------------
  cat("Fixed points, R0 and stability...\n")
  absent <- list(State_0 = c("Forage", "Predator"), State_A = character(),
                 State_B = "Forage", State_C = "Predator")
  res$states <- do.call(rbind, lapply(names(states), function(nm) {
    p <- states[[nm]]
    b <- getBiomass(p)
    r0 <- get_R0(p)
    ev_res <- if (nm == "State_0") -min(resource_rate(p)[resource_capacity(p) > 0]) else
      leading_eigenvalue(p, F0)
    inv <- if (length(absent[[nm]])) {
      max(Re(dynamic_eigenvalues(p, F0, invaders = absent[[nm]], ref = params_A)))
    }
    data.frame(
      State = nm,
      Forage_Biomass = unname(b["Forage"]), Predator_Biomass = unname(b["Predator"]),
      R0_Forage = unname(r0["Forage"]), R0_Predator = unname(r0["Predator"]),
      drift = if (nm == "State_0") 0 else biomass_drift(p, F0),
      lambda_resident = ev_res,
      lambda_full = max(ev_res, inv)
    )
  }))
  print(res$states, row.names = FALSE)

  cat("Saddle between States A and C...\n")
  saddle <- find_saddle(params_A, params_C)
  res$saddle <- list(biomass = saddle$biomass, leading = saddle$leading,
                     n_unstable = saddle$n_unstable,
                     residual = saddle$convergence$residual,
                     converged = saddle$convergence$converged)
  print(res$saddle)

  # -------------------------------------------------------------------------
  # Mechanism: larval predation versus food, diet, survival
  # -------------------------------------------------------------------------
  res$R0_decomposition <- c(
    A = R0_mix(params_A, params_A), C = R0_mix(params_C, params_C),
    growthA_mortC = R0_mix(params_A, params_C), growthC_mortA = R0_mix(params_C, params_A)
  )
  wv <- w(params_A)
  surv <- function(p) {
    k <- wv <= species_params(p)["Predator", "w_mat"]
    exp(-sum((getMort(p)["Predator", ] / getEGrowth(p)["Predator", ] * dw(p))[k]))
  }
  res$survival_to_maturity <- c(A = surv(params_A), C = surv(params_C))
  diet <- getDiet(params_A, proportion = FALSE)
  intake <- colSums(sweep(diet["Predator", , ], 1, initialN(params_A)["Predator", ] * dw(params_A), "*"))
  res$predator_diet_A <- intake / sum(intake)
  pm_A <- getPredMort(params_A)["Predator", ]
  pm_C <- getPredMort(params_C)["Predator", ]
  larval <- wv <= 1
  res$larval_predmort <- list(max_A = max(pm_A[larval]), max_C = max(pm_C[larval]),
                              w_at_max_C = wv[which.max(pm_C)])

  # -------------------------------------------------------------------------
  # Plot 1: Size Spectra and Predation Mortality across States A, B, and C
  # -------------------------------------------------------------------------
  p_spec_A <- plotSpectra(params_A, biomass = TRUE, per_log_size = TRUE, wlim = c(1e-4, 3000)) +
    ggtitle("State A: Coexistence (Cultivated)") + theme_bw(base_size = 10)
  p_spec_B <- plotSpectra(params_B, biomass = TRUE, per_log_size = TRUE, wlim = c(1e-4, 3000)) +
    ggtitle("State B: Plankton + Predator") + theme_bw(base_size = 10)
  p_spec_C <- plotSpectra(params_C, biomass = TRUE, per_log_size = TRUE, wlim = c(1e-4, 3000)) +
    ggtitle("State C: Plankton + Forage (Uncultivated)") + theme_bw(base_size = 10)

  df_pm <- rbind(
    data.frame(w = wv, PredMort = as.numeric(pm_A), State = "State A (Coexistence)"),
    data.frame(w = wv, PredMort = as.numeric(pm_C), State = "State C (Forage only)")
  )
  p_mort_comp <- ggplot(df_pm[df_pm$w <= 10, ], aes(x = w, y = PredMort, colour = State, linetype = State)) +
    geom_line(linewidth = 1) +
    scale_x_log10() +
    scale_colour_manual(values = c("State A (Coexistence)" = "#009E73", "State C (Forage only)" = "#CC79A7")) +
    labs(
      title = "Depensatory Bottleneck: Predation Mortality on Predator",
      x = "Predator weight w [g]",
      y = "Predation mortality [1/year]"
    ) +
    theme_bw(base_size = 10)

  fig1 <- (p_spec_A | p_spec_B) / (p_spec_C | p_mort_comp)
  ggsave("plots/spectra_states_ABC.png", fig1, width = 11, height = 7.5, dpi = 150)

  # -------------------------------------------------------------------------
  # Plot 2: Community Assembly / Invasion Experiments (F = 0)
  # -------------------------------------------------------------------------
  cat("Invasion experiments...\n")
  p_inv_B <- params_B
  initialN(p_inv_B)["Forage", ] <- initialN(params_A)["Forage", ] * 1e-3
  sim_inv_B <- project(p_inv_B, t_max = 150, dt = 0.1, t_save = 1, method = "tr_bdf2", progress_bar = FALSE)

  p_inv_C <- params_C
  initialN(p_inv_C)["Predator", ] <- initialN(params_A)["Predator", ] * 1e-2
  sim_inv_C <- project(p_inv_C, t_max = 150, dt = 0.1, t_save = 1, method = "tr_bdf2", progress_bar = FALSE)
  bC <- getBiomass(sim_inv_C)[, "Predator"]
  res$invasion_C_decay_rate <- unname(log(bC[151] / bC[51]) / 100)
  bB <- getBiomass(sim_inv_B)
  res$invasion_B_final <- bB[nrow(bB), ]

  p_plot_invB <- plotBiomass(sim_inv_B) +
    ggtitle("A. Forage Fish Invading State B -> Reaches State A") +
    theme_bw(base_size = 10)
  p_plot_invC <- plotBiomass(sim_inv_C) +
    ggtitle("B. Predator Invading State C -> Fails (Depensation)") +
    theme_bw(base_size = 10)

  fig2 <- p_plot_invB | p_plot_invC
  ggsave("plots/assembly_invasions.png", fig2, width = 11, height = 4.5, dpi = 150)

  # -------------------------------------------------------------------------
  # Plot 3: Hysteresis loops along F_Pred at F_Forage = 0 and 0.35
  # -------------------------------------------------------------------------
  cat("Tracing equilibrium branches along F_Pred...\n")
  Fp_seq <- seq(0, 0.34, by = 0.01)
  seed <- initialN(params_A)["Predator", ] * 1e-6
  branch_jobs <- expand.grid(Ff = c(0, 0.35), dir = c("up", "down"), stringsAsFactors = FALSE)
  branches <- parallel::mclapply(seq_len(nrow(branch_jobs)), function(i) {
    Ff <- branch_jobs$Ff[i]
    if (branch_jobs$dir[i] == "up") {
      d <- trace_branch(params_A, Ff, Fp_seq)
      d$Sweep <- "Increasing F_Pred (from State A)"
    } else {
      d <- trace_branch(params_C, Ff, rev(Fp_seq), seed = seed, t_max = 1000)
      d$Sweep <- "Decreasing F_Pred (from State C + seed)"
    }
    d
  }, mc.cores = n_cores)
  df_branch <- do.call(rbind, branches)
  res$branches <- df_branch

  df_branch_long <- rbind(
    transform(df_branch, Species = "Forage", Biomass = Forage),
    transform(df_branch, Species = "Predator", Biomass = Predator)
  )
  df_branch_long$Panel <- sprintf("F_Forage = %.2f", df_branch_long$F_Forage)
  df_branch_long$Sweep <- factor(df_branch_long$Sweep, levels = c(
    "Increasing F_Pred (from State A)", "Decreasing F_Pred (from State C + seed)"))
  fig3 <- ggplot(df_branch_long, aes(x = F_Pred, y = Biomass, colour = Species, linetype = Sweep)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.2) +
    facet_wrap(~Panel) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                       breaks = c(0, 1, 10, 100, 1000)) +
    scale_colour_manual(values = c(Forage = "#0072B2", Predator = "#D55E00")) +
    labs(
      title = "Equilibrium biomass along increasing and decreasing predator fishing mortality",
      x = "Predator fishing mortality F_Pred [1/year]",
      y = "Equilibrium biomass [g] (pseudo-log scale)",
      linetype = NULL
    ) +
    theme_bw(base_size = 10.5) +
    theme(legend.position = "bottom", legend.box = "vertical")
  ggsave("plots/predator_exploitation_bifurcation.png", fig3, width = 10.5, height = 5.8, dpi = 150)

  # -------------------------------------------------------------------------
  # Plot 4: Multi-Phase Time Series (Collapse, Irreversibility, & Forage-Harvest Recovery)
  #   Phase 1 (0-50):    Unfished State A (F_Pred = 0, F_Forage = 0)
  #   Phase 2 (50-200):  Predator overfished (F_Pred = 0.25, F_Forage = 0) -> Collapse to State C
  #   Phase 3 (200-320): Predator moratorium + re-introduction seed (F_Pred = 0, F_Forage = 0) -> Cannot re-enter
  #   Phase 4 (320-470): Harvest Forage fish + seed (F_Pred = 0, F_Forage = 0.40, above
  #                      the recovery threshold F_forage_rec) -> Predator recovers
  #   Phase 5 (470-570): All fishing stopped (F_Pred = 0, F_Forage = 0) -> Coexistence State A self-sustains
  # -------------------------------------------------------------------------
  cat("Management time series...\n")
  f_over <- 0.25
  f_forage_restore <- 0.40
  sim1 <- project(params_A, effort = F0, t_max = 50, dt = 0.1, t_save = 1, method = "tr_bdf2", progress_bar = FALSE)
  sim2 <- project(finalParams(sim1), effort = c(Forage = 0, Predator = f_over), t_max = 150, dt = 0.1, t_save = 1, method = "tr_bdf2", progress_bar = FALSE)

  p_moratorium <- finalParams(sim2)
  initialN(p_moratorium)["Predator", ] <- initialN(params_A)["Predator", ] * 5e-3
  sim3 <- project(p_moratorium, effort = F0, t_max = 120, dt = 0.1, t_save = 1, method = "tr_bdf2", progress_bar = FALSE)

  p_harvest_forage <- finalParams(sim3)
  initialN(p_harvest_forage)["Predator", ] <- initialN(params_A)["Predator", ] * 5e-3
  sim4 <- project(p_harvest_forage, effort = c(Forage = f_forage_restore, Predator = 0), t_max = 150, dt = 0.1, t_save = 1, method = "tr_bdf2", progress_bar = FALSE)

  sim5 <- project(finalParams(sim4), effort = F0, t_max = 100, dt = 0.1, t_save = 1, method = "tr_bdf2", progress_bar = FALSE)

  extract_ts <- function(sim, t_offset) {
    b <- getBiomass(sim)
    t_vec <- as.numeric(rownames(b)) + t_offset
    rbind(
      data.frame(Year = t_vec, Species = "Forage", Biomass = as.numeric(b[, "Forage"])),
      data.frame(Year = t_vec, Species = "Predator", Biomass = as.numeric(b[, "Predator"]))
    )
  }

  df_ts <- rbind(
    extract_ts(sim1, 0),
    extract_ts(sim2, 50)[-c(1, 2), ],
    extract_ts(sim3, 200)[-c(1, 2), ],
    extract_ts(sim4, 320)[-c(1, 2), ],
    extract_ts(sim5, 470)[-c(1, 2), ]
  )
  b4 <- getBiomass(sim4)
  b5 <- getBiomass(sim5)
  res$hysteresis <- list(
    f_over = f_over, f_forage_restore = f_forage_restore,
    end_phase2 = getBiomass(sim2)[nrow(getBiomass(sim2)), ],
    end_phase4 = b4[nrow(b4), ], end_phase5 = b5[nrow(b5), ],
    year_pred_above_half_A = 320 + as.numeric(rownames(b4))[
      which(b4[, "Predator"] > 0.5 * getBiomass(params_A)["Predator"])[1]]
  )

  fig4 <- ggplot(df_ts, aes(x = Year, y = Biomass, colour = Species)) +
    geom_line(linewidth = 1) +
    scale_y_log10(limits = c(1e-1, 5e3)) +
    scale_colour_manual(values = c(Forage = "#0072B2", Predator = "#D55E00")) +
    geom_vline(xintercept = c(50, 200, 320, 470), linetype = "dashed", colour = "grey40") +
    annotate("text", x = 25,  y = 3800, label = "1. Unfished\n(State A)", size = 3) +
    annotate("text", x = 125, y = 3800, label = sprintf("2. Predator fished\n(F_Pred = %.2f)\n-> Collapse to C", f_over), size = 3) +
    annotate("text", x = 260, y = 3800, label = "3. Moratorium + Seed\n(F_Pred = 0, F_For = 0)\n-> Cannot re-enter", size = 3) +
    annotate("text", x = 395, y = 3800, label = sprintf("4. Harvest Forage + Seed\n(F_Forage = %.2f)\n-> Predator recovers", f_forage_restore), size = 3) +
    annotate("text", x = 520, y = 3800, label = "5. Unfished\n(State A\nrestored)", size = 3) +
    labs(
      title = "Cultivation-Depensation Hysteresis: Collapse, Failed Moratorium, and Recovery via Forage Harvesting",
      x = "Time [years]",
      y = "Biomass [g] (log10 scale)"
    ) +
    theme_bw(base_size = 10.5)
  ggsave("plots/hysteresis_and_recovery.png", fig4, width = 11, height = 5.5, dpi = 150)

  # -------------------------------------------------------------------------
  # Plot 5: Predator Invasion Fitness in State C & 2D Exploitation Regime Diagram
  # -------------------------------------------------------------------------
  cat("Predator R0 in State C along F_Forage...\n")
  Ff_fine <- seq(0, 0.70, by = 0.025)
  r0_pred_C <- numeric(length(Ff_fine))
  b_for_C <- numeric(length(Ff_fine))
  p_C_curr <- params_C
  for (i in seq_along(Ff_fine)) {
    eff <- c(Forage = Ff_fine[i], Predator = 0)
    p_C_curr <- settle(p_C_curr, eff)
    b_for_C[i] <- getBiomass(p_C_curr)["Forage"]
    r0_pred_C[i] <- get_R0(p_C_curr, effort = eff, species = "Predator")["Predator"]
  }
  df_inv_curve <- data.frame(F_Forage = Ff_fine, Forage_Biomass_C = b_for_C, R0_Pred_C = r0_pred_C)
  res$r0_curve <- df_inv_curve
  f_r0 <- function(Ff) {
    p <- settle(params_C, c(Forage = Ff, Predator = 0))
    get_R0(p, effort = c(Forage = Ff, Predator = 0), species = "Predator")["Predator"] - 1
  }
  k <- which(diff(sign(r0_pred_C - 1)) != 0)[1]
  res$F_forage_inv <- uniroot(f_r0, Ff_fine[c(k, k + 1)], tol = 1e-5)$root
  res$B_forage_C_at_inv <- unname(getBiomass(settle(params_C, c(Forage = res$F_forage_inv, Predator = 0)))["Forage"])

  cat("Mapping regime boundaries in the (F_Forage, F_Pred) plane...\n")
  Ff_map <- seq(0, 0.70, by = 0.025)
  bounds <- do.call(rbind, parallel::mclapply(Ff_map, function(Ff) {
    map_regime(params_A, params_C, Ff)
  }, mc.cores = n_cores))
  res$bounds <- bounds
  print(bounds, row.names = FALSE)
  res$F_crit <- bounds$F_up[bounds$F_Forage == 0]
  hyst <- bounds$F_up - bounds$F_low
  # Where the hysteresis vanishes (a cusp, at which the two folds merge)
  res$cusp_F_Forage <- if (any(hyst <= 0.003)) min(bounds$F_Forage[hyst <= 0.003]) else NA
  # Where the invasion boundary changes from subcritical to supercritical
  res$tcsn_F_Forage <- if (any(bounds$A_low)) min(bounds$F_Forage[bounds$A_low]) else NA

  res$F_forage_rec <- find_F_forage_rec(params_A, params_C, lower = res$F_forage_inv)

  # The boundaries vary smoothly with F_Forage, so interpolate them onto a finer
  # grid for drawing; F_inv starts from zero at F_Forage,inv.
  Ff_draw <- seq(0, 0.7, by = 0.0025)
  has_inv <- !is.na(bounds$F_inv)
  inv_x <- c(res$F_forage_inv, bounds$F_Forage[has_inv])
  inv_y <- c(0, bounds$F_inv[has_inv])
  bounds_draw <- data.frame(
    F_Forage = Ff_draw,
    F_up = approx(bounds$F_Forage, bounds$F_up, xout = Ff_draw)$y,
    F_inv = ifelse(Ff_draw < res$F_forage_inv, NA, approx(inv_x, inv_y, xout = Ff_draw, rule = 2)$y),
    A_low = approx(bounds$F_Forage, as.numeric(bounds$A_low), xout = Ff_draw, method = "constant", f = 1)$y > 0.5
  )
  low_x <- c(bounds$F_Forage, res$F_forage_inv)
  low_y <- c(ifelse(bounds$A_low | is.na(bounds$F_inv), bounds$F_low, NA), 0)
  ok <- !is.na(low_y)
  bounds_draw$F_low <- ifelse(bounds_draw$A_low,
                              approx(low_x[ok], low_y[ok], xout = Ff_draw, rule = 2)$y,
                              pmax(0, ifelse(is.na(bounds_draw$F_inv), 0, bounds_draw$F_inv)))
  regimes <- classify_regimes(bounds_draw, seq(0, 0.6, by = 0.0025))
  regime_cols <- c(
    "Coexistence only" = "#009E73",
    "Bistable: A and C" = "#E69F00",
    "Bistable: A and low-predator state" = "#F0E442",
    "Low-predator coexistence only" = "#CC79A7",
    "State C only" = "#56B4E9"
  )
  p_r0 <- ggplot(df_inv_curve, aes(x = F_Forage, y = R0_Pred_C)) +
    geom_line(linewidth = 1.1, colour = "#D55E00") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
    geom_vline(xintercept = res$F_forage_inv, linetype = "dotted", colour = "#0072B2") +
    annotate("text", x = res$F_forage_inv + 0.007, y = 0.65,
             label = sprintf("Invasion threshold\nF_Forage = %.3f /yr\n(R0_Pred = 1)", res$F_forage_inv),
             hjust = 0, size = 3.1, colour = "#0072B2") +
    labs(
      title = "A. Predator Invasion Fitness in State C",
      x = "Forage fishing mortality F_Forage [1/year]",
      y = "Predator invasion R0 in State C"
    ) +
    theme_bw(base_size = 10.5)

  cusp_draw <- ifelse(is.na(res$cusp_F_Forage), Inf, res$cusp_F_Forage)
  boundary_types <- c("Upper limit of abundant coexistence (fold)" = "solid",
                      "Lower limit of low-predator state" = "dotted",
                      "Predator invasion of State C (R0_Pred = 1)" = "dashed")
  bound_lines <- rbind(
    data.frame(F_Forage = bounds_draw$F_Forage,
               F_Pred = ifelse(bounds_draw$F_Forage <= cusp_draw, bounds_draw$F_up, NA),
               Boundary = names(boundary_types)[1]),
    data.frame(F_Forage = bounds_draw$F_Forage,
               F_Pred = ifelse(bounds_draw$A_low & bounds_draw$F_Forage <= cusp_draw, bounds_draw$F_low, NA),
               Boundary = names(boundary_types)[2]),
    data.frame(F_Forage = bounds_draw$F_Forage, F_Pred = bounds_draw$F_inv,
               Boundary = names(boundary_types)[3])
  )
  bound_lines <- bound_lines[!is.na(bound_lines$F_Pred) & bound_lines$F_Pred > 0, ]
  bound_lines$Boundary <- factor(bound_lines$Boundary, levels = names(boundary_types))
  p_phase <- ggplot() +
    geom_tile(data = regimes, aes(x = F_Forage, y = F_Pred, fill = Regime),
              width = 0.0025, height = 0.0025, alpha = 0.55) +
    geom_line(data = bound_lines, aes(x = F_Forage, y = F_Pred, linetype = Boundary), linewidth = 0.8) +
    scale_fill_manual(values = regime_cols) +
    scale_linetype_manual(values = boundary_types) +
    coord_cartesian(xlim = c(0, 0.70), ylim = c(0, 0.6), expand = FALSE) +
    labs(
      title = "B. Two-Parameter Exploitation Regime Diagram",
      x = "Forage fishing mortality F_Forage [1/year]",
      y = "Predator fishing mortality F_Pred [1/year]",
      fill = NULL, linetype = NULL
    ) +
    guides(fill = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
    theme_bw(base_size = 10.5) +
    theme(legend.position = "right")

  fig5 <- (p_r0 | p_phase) + plot_layout(widths = c(1, 1.6))
  ggsave("plots/regime_diagram_2D.png", fig5, width = 14, height = 5.6, dpi = 150)

  # -------------------------------------------------------------------------
  # Robustness: grid resolution, resource cutoff, intraspecific predation
  # -------------------------------------------------------------------------
  cat("Robustness checks...\n")
  run_jobs <- function(jobs, thresholds) {
    cols <- c("B_For_A", "B_Pred_A", "B_For_C", "R0_For_B", "R0_Pred_C", "Forage_in_diet",
              if (thresholds) c("F_crit", "F_forage_inv", "F_forage_rec"))
    do.call(rbind, parallel::mclapply(jobs, function(j) {
      tryCatch(do.call(variant_summary, c(j, thresholds = thresholds)), error = function(e) {
        message("Variant ", j[[1]], " failed: ", conditionMessage(e))
        cbind(data.frame(Variant = j[[1]]), as.data.frame(as.list(setNames(rep(NA_real_, length(cols)), cols))))
      })
    }, mc.cores = n_cores))
  }
  # Numerical resolution and model structure, with the fishing thresholds
  variants <- list(
    list("Reference (320 size classes)"),
    list("160 size classes", no_w = 160),
    list("640 size classes", no_w = 640),
    list("Original model (first-order, 80 classes, no diffusion)",
         no_w = 80, flux = "upwind", c_diff = 0, min_w_pp = 1e-6,
         selectivity = "knife_edge", ogive_l25 = NA, taper = NA,
         w_max_factor = 1, z_sen = 0),
    list("Resource tapered at 1 g", w_pp_cutoff = 1)
  )
  res$robustness <- run_jobs(variants, thresholds = TRUE)
  print(res$robustness, row.names = FALSE)

  # Strength of the physical size diffusion. The centred flux stays monotone
  # only while the cell Peclet number dx / c_diff is below about 1; for weaker
  # diffusion the van Leer limiter keeps the solution non-negative.
  c_grid <- c(0.02, 0.05, 0.1, 0.2, 0.3)
  dx_ref <- log(w(params_A)[2] / w(params_A)[1])
  diff_jobs <- lapply(c_grid, function(cd) {
    list(sprintf("c = %.2f", cd), c_diff = cd,
         flux = if (dx_ref / cd < 1) "centred" else "van_leer")
  })
  res$diffusion_scan <- cbind(c_diff = c_grid, run_jobs(diff_jobs, thresholds = TRUE))
  print(res$diffusion_scan, row.names = FALSE)

  # Strength of cannibalism
  theta_grid <- c(0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 0.5, 0.95)
  scan_jobs <- function(param) {
    lapply(theta_grid, function(th) {
      j <- list(sprintf("%s = %.3f", param, th))
      j[[param]] <- th
      j
    })
  }
  rob <- run_jobs(c(scan_jobs("theta_FF"), scan_jobs("theta_PP")), thresholds = FALSE)
  ng <- length(theta_grid)
  res$theta_scan <- rbind(
    cbind(Parameter = "theta_FF", theta = theta_grid, rob[seq_len(ng), ]),
    cbind(Parameter = "theta_PP", theta = theta_grid, rob[ng + seq_len(ng), ])
  )
  # Interaction strength at which the predator first becomes able to invade State C
  theta_crit <- function(param) {
    d <- res$theta_scan[res$theta_scan$Parameter == param, ]
    k <- which(diff(sign(log(d$R0_Pred_C))) != 0)[1]
    if (is.na(k)) return(NA_real_)
    approx(log(d$R0_Pred_C[k:(k + 1)]), d$theta[k:(k + 1)], xout = 0)$y
  }
  res$theta_FF_crit <- theta_crit("theta_FF")
  res$theta_PP_crit <- theta_crit("theta_PP")

  ts <- res$theta_scan
  ts$Parameter <- factor(ts$Parameter, levels = c("theta_FF", "theta_PP"),
                         labels = c("Forage cannibalism theta_For,For", "Predator cannibalism theta_Pred,Pred"))
  crit <- data.frame(Parameter = levels(ts$Parameter), theta = c(res$theta_FF_crit, res$theta_PP_crit))
  fig6 <- ggplot(ts, aes(x = theta, y = R0_Pred_C)) +
    geom_line(linewidth = 1, colour = "#D55E00") +
    geom_point(size = 2, colour = "#D55E00") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
    geom_vline(data = crit[!is.na(crit$theta), ], aes(xintercept = theta),
               linetype = "dotted", colour = "#0072B2") +
    facet_wrap(~Parameter) +
    scale_y_log10() +
    labs(
      title = "Predator invasion R0 in State C versus the strength of cannibalism",
      x = "Interaction strength (the other entries of the interaction matrix unchanged)",
      y = "Predator invasion R0 in State C (log scale)"
    ) +
    theme_bw(base_size = 10.5)
  ggsave("plots/robustness_cannibalism.png", fig6, width = 10, height = 4.5, dpi = 150)

  ds <- res$diffusion_scan
  ds_long <- rbind(
    data.frame(c_diff = ds$c_diff, Value = ds$F_crit, Quantity = "F_crit (collapse, F_Forage = 0)"),
    data.frame(c_diff = ds$c_diff, Value = ds$F_forage_inv, Quantity = "F_Forage,inv (predator can invade C)"),
    data.frame(c_diff = ds$c_diff, Value = ds$F_forage_rec, Quantity = "F_Forage,rec (predator recovers fully)")
  )
  p_thr <- ggplot(ds_long, aes(x = c_diff, y = Value, colour = Quantity)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_log10() +
    labs(title = "A. Fishing thresholds", x = "Size diffusion c = d sigma^2 / d ln w",
         y = "Fishing mortality [1/year]", colour = NULL) +
    theme_bw(base_size = 10.5) +
    theme(legend.position = "bottom", legend.direction = "vertical")
  p_r0d <- ggplot(ds, aes(x = c_diff, y = R0_Pred_C)) +
    geom_line(linewidth = 1, colour = "#D55E00") +
    geom_point(size = 2, colour = "#D55E00") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
    scale_x_log10() +
    labs(title = "B. Predator invasion R0 in State C", x = "Size diffusion c = d sigma^2 / d ln w",
         y = "R0_Pred(State C)") +
    theme_bw(base_size = 10.5)
  ggsave("plots/robustness_diffusion.png", p_thr | p_r0d, width = 10, height = 4.8, dpi = 150)

  saveRDS(res, "results/cultivation_results.rds")
  cat(sprintf("\nPredator collapse threshold (F_Forage = 0): F_crit = %.3f /year\n", res$F_crit))
  cat("All plots saved to plots/, results to results/cultivation_results.rds,\n")
  cat("and model saved to cultivation_depensation_params.rds.\n")
}
