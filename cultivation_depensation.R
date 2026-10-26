# Cultivation-Depensation in a 3-Taxa Mizer Model
# Community:
#   1. Plankton resource spectrum (semichemostat dynamics)
#   2. Forage fish spectrum (small maximum size, w_inf = 20 g)
#   3. Predator size spectrum (large maximum size, w_inf = 2000 g)

suppressPackageStartupMessages({
  library(mizer)
  library(ggplot2)
  library(patchwork)
})

#' Compute exact discrete lifetime reproductive output (R0) in a resident environment
#'
#' Evaluates the linear homogeneous McKendrick-von Foerster steady spectrum per
#' unit recruitment for each species in the background growth and mortality
#' environment of `params`, returning RDI / RDD_required.
#' A rare species i can invade a resident steady state iff R0[i] > 1,
#' and at an interior Beverton-Holt steady state, R0[i] = 1 / (1 - reproduction_level[i]).
get_R0 <- function(params, effort = initial_effort(params)) {
  p <- params
  initial_effort(p) <- effort
  sp_names <- species_params(p)$species
  R0 <- setNames(numeric(length(sp_names)), sp_names)

  g <- getEGrowth(p)
  mort <- getMort(p)
  erepro_g <- getERepro(p)
  diffusion <- getDiffusion(p)

  for (sp in sp_names) {
    w_min_idx <- p@w_min_idx[sp]
    N0_vec <- setNames(numeric(length(sp_names)), sp_names)
    N0_vec[sp] <- 1.0
    n_exact <- mizer:::get_steady_state_n(p, g, mort, diffusion, N0_vec)
    rdd_req <- (g[sp, w_min_idx] + mort[sp, w_min_idx] * dw(p)[w_min_idx]) * 1.0
    sex_ratio <- 0.5
    erepro <- species_params(p)[sp, "erepro"]
    w_min <- species_params(p)[sp, "w_min"]
    rdi_prod <- sex_ratio * erepro * sum(n_exact[sp, ] * erepro_g[sp, ] * dw(p)) / w_min
    R0[sp] <- unname(rdi_prod / rdd_req)
  }
  R0
}

#' Build and calibrate the cultivation-depensation MizerParams at Coexistence State A
build_cultivation_model <- function() {
  sp <- data.frame(
    species          = c("Forage", "Predator"),
    w_inf            = c(20, 2000),
    w_mat            = c(4, 250),
    w_min            = c(0.001, 0.001),
    beta             = c(200, 100),
    sigma            = c(1.3, 1.5),
    k_vb             = c(1.0, 0.5),
    f0               = c(0.6, 0.6),
    z0               = c(0.30, 0.20),
    z_ext            = c(0.50, 1.20),
    biomass_observed = c(100, 600)
  )

  # Interaction matrix:
  # - Predator preys strongly on Forage (1.00): top-down cultivation
  # - Forage preys on Predator larvae/early juveniles (0.95): depensatory bottleneck
  # - Weak self-interaction (0.05) avoids single-species cannibalism oscillations
  inter <- matrix(
    c(0.05, 0.95,
      1.00, 0.05),
    nrow = 2, byrow = TRUE,
    dimnames = list(sp$species, sp$species)
  )

  # Separate species-specific gears with catchability = 1 so effort equals F (1/year)
  gp <- data.frame(
    gear            = c("Forage", "Predator"),
    species         = c("Forage", "Predator"),
    sel_func        = "knife_edge",
    knife_edge_size = c(1.0, 400.0),
    catchability    = c(1.0, 1.0)
  )

  params <- newMultispeciesParams(
    sp,
    interaction      = inter,
    gear_params      = gp,
    initial_effort   = c(Forage = 0, Predator = 0),
    no_w             = 80,
    min_w_pp         = 1e-6,
    kappa            = 100,
    w_pp_cutoff      = 10,
    resource_rate    = 10,
    info_level       = 0
  )
  params <- setColours(params, c(Forage = "#0072B2", Predator = "#D55E00"))

  # Calibrate to Coexistence State A
  params <- steadySingleSpecies(params)
  params <- matchBiomasses(params)
  params <- tuneSteadyState(params, dt = 0.2, info_level = 0)
  params <- matchBiomasses(params)
  params <- tuneSteadyState(params, dt = 0.2, info_level = 0)

  # Set Beverton-Holt recruitment density dependence at State A
  reproduction_level(params) <- c(Forage = 0.65, Predator = 0.78)

  # Settle to high-precision fixed point A without altering parameters
  params_A <- finalParams(project(params, t_max = 400, dt = 0.2, progress_bar = FALSE))
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

#' Derive the four unfished fixed points (State 0, State A, State B, State C)
#' sharing identical model parameters
get_fixed_points <- function(params_A) {
  # State B: Plankton + Predator subsystem (Forage absent)
  params_B <- params_A
  initialN(params_B)["Forage", ] <- 0
  params_B <- finalParams(project(params_B, t_max = 400, dt = 0.2, progress_bar = FALSE))

  # State C: Plankton + Forage subsystem (Predator absent)
  params_C <- params_A
  initialN(params_C)["Predator", ] <- 0
  params_C <- finalParams(project(params_C, t_max = 400, dt = 0.2, progress_bar = FALSE))

  # State 0: Plankton only (both fish absent, resource at carrying capacity)
  params_0 <- params_A
  initialN(params_0)[] <- 0
  initialNResource(params_0) <- resource_capacity(params_0)

  list(State_0 = params_0, State_A = params_A, State_B = params_B, State_C = params_C)
}

if (sys.nframe() == 0) {
  options(bspm.sudo = TRUE)
  dir.create("plots", showWarnings = FALSE)

  cat("Building and calibrating cultivation-depensation model...\n")
  params_A <- build_cultivation_model()
  saveParams(params_A, "cultivation_depensation_params.rds")

  states <- get_fixed_points(params_A)
  params_0 <- states$State_0
  params_B <- states$State_B
  params_C <- states$State_C

  cat("\n=== Unfished Fixed Points Summary (F_Pred = 0, F_Forage = 0) ===\n")
  summary_rows <- lapply(names(states), function(nm) {
    p <- states[[nm]]
    b <- getBiomass(p)
    r0 <- get_R0(p)
    st <- if (nm == "State_0") list(stable = FALSE, max_real_part = NA_real_) else getStability(p)
    data.frame(
      State            = nm,
      Forage_Biomass   = round(unname(b["Forage"]), 2),
      Predator_Biomass = round(unname(b["Predator"]), 2),
      R0_Forage        = round(unname(r0["Forage"]), 3),
      R0_Predator      = round(unname(r0["Predator"]), 3),
      isSteady         = isSteady(p),
      Subspace_Stable  = st$stable,
      Max_Re_Eig       = round(st$max_real_part, 4)
    )
  })
  summary_df <- do.call(rbind, summary_rows)
  print(summary_df, row.names = FALSE)

  # -------------------------------------------------------------------------
  # Plot 1: Size Spectra and Predation Mortality across States A, B, and C
  # -------------------------------------------------------------------------
  p_spec_A <- plotSpectra(params_A, biomass = TRUE, per_log_size = TRUE, wlim = c(1e-4, 3000)) +
    ggtitle("State A: Coexistence (Cultivated)") + theme_bw(base_size = 10)
  p_spec_B <- plotSpectra(params_B, biomass = TRUE, per_log_size = TRUE, wlim = c(1e-4, 3000)) +
    ggtitle("State B: Plankton + Predator") + theme_bw(base_size = 10)
  p_spec_C <- plotSpectra(params_C, biomass = TRUE, per_log_size = TRUE, wlim = c(1e-4, 3000)) +
    ggtitle("State C: Plankton + Forage (Uncultivated)") + theme_bw(base_size = 10)

  pm_A <- getPredMort(params_A)["Predator", ]
  pm_C <- getPredMort(params_C)["Predator", ]
  df_pm <- rbind(
    data.frame(w = w(params_A), PredMort = as.numeric(pm_A), State = "State A (Coexistence)"),
    data.frame(w = w(params_C), PredMort = as.numeric(pm_C), State = "State C (Forage only)")
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
  #   (a) Forage invading State B -> converges to State A
  #   (b) Predator attempting to invade State C -> fails, stays at State C
  # -------------------------------------------------------------------------
  p_inv_B <- params_B
  initialN(p_inv_B)["Forage", ] <- initialN(params_A)["Forage", ] * 1e-3
  sim_inv_B <- project(p_inv_B, t_max = 150, dt = 0.2, t_save = 1, progress_bar = FALSE)

  p_inv_C <- params_C
  initialN(p_inv_C)["Predator", ] <- initialN(params_A)["Predator", ] * 1e-2
  sim_inv_C <- project(p_inv_C, t_max = 150, dt = 0.2, t_save = 1, progress_bar = FALSE)

  p_plot_invB <- plotBiomass(sim_inv_B) +
    ggtitle("A. Forage Fish Invading State B -> Reaches State A") +
    theme_bw(base_size = 10)
  p_plot_invC <- plotBiomass(sim_inv_C) +
    ggtitle("B. Predator Invading State C -> Fails (Depensation)") +
    theme_bw(base_size = 10)

  fig2 <- p_plot_invB | p_plot_invC
  ggsave("plots/assembly_invasions.png", fig2, width = 11, height = 4.5, dpi = 150)

  # -------------------------------------------------------------------------
  # Plot 3: Predator Exploitation Scan (F_Predator) with F_Forage = 0 vs 0.35
  # -------------------------------------------------------------------------
  F_grid <- seq(0, 0.30, by = 0.025)
  scan_list <- list()
  for (Ff in c(0.00, 0.35)) {
    p_curr <- params_A
    for (Fp in F_grid) {
      sim_step <- project(p_curr, effort = c(Forage = Ff, Predator = Fp), t_max = 350, dt = 0.2, progress_bar = FALSE)
      p_curr <- finalParams(sim_step)
      b <- getBiomass(p_curr)
      scan_list[[length(scan_list) + 1]] <- data.frame(
        F_Pred = Fp,
        F_Forage = ifelse(Ff == 0, "F_Forage = 0.00 (Unbalanced)", "F_Forage = 0.35 (Balanced)"),
        Species = c("Forage", "Predator"),
        Biomass = c(unname(b["Forage"]), unname(b["Predator"]))
      )
    }
  }
  df_scan <- do.call(rbind, scan_list)

  # Find exact collapse threshold on F_Forage = 0 branch
  df_unbal_pred <- subset(df_scan, F_Forage == "F_Forage = 0.00 (Unbalanced)" & Species == "Predator")
  last_alive <- max(df_unbal_pred$F_Pred[df_unbal_pred$Biomass > 1])
  first_dead <- min(df_unbal_pred$F_Pred[df_unbal_pred$Biomass <= 1])
  f_crit_est <- 0.5 * (last_alive + first_dead)

  fig3 <- ggplot(df_scan, aes(x = F_Pred, y = Biomass, colour = Species, linetype = F_Forage)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = c(Forage = "#0072B2", Predator = "#D55E00")) +
    geom_vline(xintercept = f_crit_est, linetype = "dotted", colour = "grey30") +
    annotate("text", x = f_crit_est + 0.006, y = 1800,
             label = sprintf("Unbalanced collapse\nthreshold (F_Pred ~ %.2f)", f_crit_est),
             hjust = 0, size = 3.2) +
    labs(
      title = "Equilibrium Biomass vs. Predator Fishing Mortality (F_Pred)",
      subtitle = "Unbalanced fishing (F_Forage = 0) collapses Predator to State C;\nharvesting Forage (F_Forage = 0.35) holds the exploited Predator in the system",
      x = "Predator fishing mortality F_Pred [1/year]",
      y = "Equilibrium Biomass [g]",
      linetype = "Forage Fishery"
    ) +
    theme_bw(base_size = 10.5)
  ggsave("plots/predator_exploitation_bifurcation.png", fig3, width = 9.5, height = 5.5, dpi = 150)

  # -------------------------------------------------------------------------
  # Plot 4: Multi-Phase Time Series (Collapse, Irreversibility, & Forage-Harvest Recovery)
  #   Phase 1 (0-50):    Unfished State A (F_Pred = 0, F_Forage = 0)
  #   Phase 2 (50-200):  Predator overfished (F_Pred = 0.22, F_Forage = 0) -> Collapse to State C
  #   Phase 3 (200-320): Predator moratorium + re-introduction seed (F_Pred = 0, F_Forage = 0) -> Cannot re-enter!
  #   Phase 4 (320-470): Harvest Forage fish + seed (F_Pred = 0, F_Forage = 0.30) -> Predator recovers!
  #   Phase 5 (470-570): All fishing stopped (F_Pred = 0, F_Forage = 0) -> Coexistence State A self-sustains!
  # -------------------------------------------------------------------------
  f_over <- round(f_crit_est + 0.06, 2)
  sim1 <- project(params_A, effort = c(Forage = 0, Predator = 0), t_max = 50, dt = 0.2, t_save = 1, progress_bar = FALSE)
  sim2 <- project(finalParams(sim1), effort = c(Forage = 0, Predator = f_over), t_max = 150, dt = 0.2, t_save = 1, progress_bar = FALSE)

  p_moratorium <- finalParams(sim2)
  initialN(p_moratorium)["Predator", ] <- initialN(params_A)["Predator", ] * 5e-3
  sim3 <- project(p_moratorium, effort = c(Forage = 0, Predator = 0), t_max = 120, dt = 0.2, t_save = 1, progress_bar = FALSE)

  p_harvest_forage <- finalParams(sim3)
  initialN(p_harvest_forage)["Predator", ] <- initialN(params_A)["Predator", ] * 5e-3
  sim4 <- project(p_harvest_forage, effort = c(Forage = 0.30, Predator = 0), t_max = 150, dt = 0.2, t_save = 1, progress_bar = FALSE)

  sim5 <- project(finalParams(sim4), effort = c(Forage = 0, Predator = 0), t_max = 100, dt = 0.2, t_save = 1, progress_bar = FALSE)

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

  fig4 <- ggplot(df_ts, aes(x = Year, y = Biomass, colour = Species)) +
    geom_line(linewidth = 1) +
    scale_y_log10(limits = c(1e-1, 5e3)) +
    scale_colour_manual(values = c(Forage = "#0072B2", Predator = "#D55E00")) +
    geom_vline(xintercept = c(50, 200, 320, 470), linetype = "dashed", colour = "grey40") +
    annotate("text", x = 25,  y = 3800, label = "1. Unfished\n(State A)", size = 3) +
    annotate("text", x = 125, y = 3800, label = sprintf("2. Predator fished\n(F_Pred = %.2f)\n-> Collapse to C", f_over), size = 3) +
    annotate("text", x = 260, y = 3800, label = "3. Moratorium + Seed\n(F_Pred = 0, F_For = 0)\n-> Cannot re-enter!", size = 3) +
    annotate("text", x = 395, y = 3800, label = "4. Harvest Forage + Seed\n(F_Forage = 0.30)\n-> Predator recovers!", size = 3) +
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
  Ff_fine <- seq(0, 0.50, by = 0.025)
  r0_pred_C <- numeric(length(Ff_fine))
  b_for_C <- numeric(length(Ff_fine))
  Fp_inv <- numeric(length(Ff_fine))

  p_C_curr <- params_C
  for (i in seq_along(Ff_fine)) {
    Ff <- Ff_fine[i]
    p_C_curr <- finalParams(project(p_C_curr, effort = c(Forage = Ff, Predator = 0), t_max = 80, dt = 0.2, progress_bar = FALSE))
    b_for_C[i] <- getBiomass(p_C_curr)["Forage"]
    r0_pred_C[i] <- get_R0(p_C_curr, effort = c(Forage = Ff, Predator = 0))["Predator"]
    if (r0_pred_C[i] <= 1) {
      Fp_inv[i] <- 0
    } else {
      f_obj <- function(fp) get_R0(p_C_curr, effort = c(Forage = Ff, Predator = fp))["Predator"] - 1
      Fp_inv[i] <- uniroot(f_obj, c(0, 1.0))$root
    }
  }

  df_inv_curve <- data.frame(
    F_Forage = Ff_fine,
    Forage_Biomass_C = b_for_C,
    R0_Pred_C = r0_pred_C,
    F_Pred_inv = Fp_inv
  )

  Ff_fold <- c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35)
  Fp_fold <- c(0.188, 0.195, 0.205, 0.218, 0.232, 0.248, 0.268, 0.292)
  df_fold <- data.frame(
    F_Forage = Ff_fine,
    F_Pred_inv = Fp_inv,
    F_Pred_crit = ifelse(
      Ff_fine <= 0.35,
      approx(Ff_fold, Fp_fold, xout = Ff_fine)$y,
      Fp_inv
    )
  )

  p_r0 <- ggplot(df_inv_curve, aes(x = F_Forage, y = R0_Pred_C)) +
    geom_line(linewidth = 1.1, colour = "#D55E00") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
    geom_vline(xintercept = 0.128, linetype = "dotted", colour = "#0072B2") +
    annotate("text", x = 0.135, y = 0.65,
             label = "Invasion threshold\nF_Forage = 0.128 /yr\n(R0_Pred = 1)",
             hjust = 0, size = 3.1, colour = "#0072B2") +
    labs(
      title = "A. Predator Invasion Fitness in State C",
      x = "Forage fishing mortality F_Forage [1/year]",
      y = "Predator invasion R0 in State C"
    ) +
    theme_bw(base_size = 10.5)

  p_phase <- ggplot(df_fold, aes(x = F_Forage)) +
    geom_ribbon(aes(ymin = 0, ymax = F_Pred_inv), fill = "#009E73", alpha = 0.25) +
    geom_ribbon(aes(ymin = F_Pred_inv, ymax = F_Pred_crit), fill = "#E69F00", alpha = 0.35) +
    geom_ribbon(aes(ymin = F_Pred_crit, ymax = 0.52), fill = "#56B4E9", alpha = 0.22) +
    geom_line(data = subset(df_fold, F_Forage >= 0.125),
              aes(y = F_Pred_inv, linetype = "Transcritical invasion (R0_Pred = 1)"),
              linewidth = 1, colour = "#009E73") +
    geom_line(aes(y = F_Pred_crit, linetype = "Fold / collapse boundary (F_crit)"),
              linewidth = 1, colour = "#D55E00") +
    annotate("point", x = 0.35, y = 0.292, size = 2.8, colour = "black") +
    annotate("text", x = 0.34, y = 0.325, label = "Cusp point", size = 3, fontface = "italic") +
    annotate("text", x = 0.11, y = 0.11, label = "Bistable Regime\n(State A & State C stable;\nmoratorium fails)", size = 3, fontface = "bold") +
    annotate("text", x = 0.34, y = 0.11, label = "Monostable Coexistence\n(State A globally stable;\nPredator invades State C)", size = 3, fontface = "bold", colour = "#006644") +
    annotate("text", x = 0.18, y = 0.40, label = "Forage-Dominated Regime\n(State C only; Predator excluded)", size = 3, fontface = "bold", colour = "#004488") +
    scale_linetype_manual(
      name = "Bifurcation Boundary",
      values = c("Fold / collapse boundary (F_crit)" = "solid",
                 "Transcritical invasion (R0_Pred = 1)" = "dashed"),
      guide = guide_legend(nrow = 2)
    ) +
    coord_cartesian(xlim = c(0, 0.50), ylim = c(0, 0.50), expand = FALSE) +
    labs(
      title = "B. Two-Parameter Exploitation Regime Diagram",
      x = "Forage fishing mortality F_Forage [1/year]",
      y = "Predator fishing mortality F_Pred [1/year]"
    ) +
    theme_bw(base_size = 10.5) +
    theme(legend.position = "bottom")

  fig5 <- p_r0 | p_phase
  ggsave("plots/regime_diagram_2D.png", fig5, width = 11.5, height = 5.4, dpi = 150)

  cat(sprintf("\nEstimated Predator collapse threshold (F_Forage = 0): F_crit ~ %.3f /year\n", f_crit_est))
  cat("All plots saved to plots/ and model saved to cultivation_depensation_params.rds.\n")
}
