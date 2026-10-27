# Export the cultivation-depensation mizer model for the browser simulator.
#
# Writes
#   interactive/model_data.js             -- window.MODEL: all rate arrays the
#                                            JavaScript engine needs, plus the
#                                            four unfished states 0, A, B, C
#   interactive/validation/reference.json -- mizer trajectories and R0 values
#                                            that validate.js compares against
#
# Run from the project root:  Rscript interactive/export_model.R

suppressPackageStartupMessages({
  library(mizer)
  library(jsonlite)
})
source("cultivation_depensation.R")  # get_fixed_points(), get_R0()

params <- readParams("cultivation_depensation_params.rds")

# The JavaScript engine re-implements mizer's default Euler/upwind scheme with
# the default rate functions. Refuse to export a model that uses anything else.
stopifnot(
  all(unlist(params@rates_funcs) ==
        c(Rates = "mizerRates", Encounter = "mizerEncounter",
          FeedingLevel = "mizerFeedingLevel",
          EReproAndGrowth = "mizerEReproAndGrowth",
          PredRate = "mizerPredRate", PredMort = "mizerPredMort",
          FMort = "mizerFMort", Mort = "mizerMort", ERepro = "mizerERepro",
          EGrowth = "mizerEGrowth", Diffusion = "mizerDiffusion",
          ResourceMort = "mizerResourceMort", RDI = "mizerRDI",
          RDD = "BevertonHoltRDD")[names(params@rates_funcs)]),
  params@resource_dynamics == "resource_semichemostat",
  second_order_w(params)$flux == "upwind",
  !isTRUE(second_order_w(params)$bin_average),
  max(getDiffusion(params)) == 0,
  max(params@ext_encounter) == 0,
  length(params@other_dynamics) == 0
)

sp <- species_params(params)
no_sp <- nrow(sp)
no_w <- length(w(params))
no_w_full <- length(w_full(params))
idx_sp <- (no_w_full - no_w + 1):no_w_full

offset0 <- idx_sp[1] - 1          # w[j] == w_full[j + offset0]

# mizer evaluates encounter and predation rate as FFT products of length
# no_w_full, i.e. as *circular* convolutions with the real-space kernels below.
# The JavaScript engine evaluates the same circular sums directly, wrap-around
# included, so that it reproduces mizer bit for bit rather than an idealised
# linear convolution. Check that here before exporting.
phi_e <- Re(t(mvfft(t(params@ft_pred_kernel_e), inverse = TRUE))) / no_w_full
phi_p <- Re(t(mvfft(t(params@ft_pred_kernel_p), inverse = TRUE))) / no_w_full
circular <- function(phi, x, i, out_idx) {
  vapply(out_idx, function(J) {
    m <- (J - seq_len(no_w_full)) %% no_w_full
    sum(phi[i, m + 1] * x)
  }, numeric(1))
}
n <- initialN(params)
n_pp <- initialNResource(params)
prey <- outer(sp$interaction_resource, n_pp)
prey[, idx_sp] <- prey[, idx_sp] + interaction_matrix(params) %*% n
prey <- sweep(prey, 2, w_full(params) * dw_full(params), "*")
enc_direct <- t(sapply(seq_len(no_sp), function(i) circular(phi_e, prey[i, ], i, idx_sp)))
enc_direct[enc_direct < 1e-18] <- 0
enc_direct <- enc_direct * params@search_vol
enc_mizer <- getEncounter(params)
rel_err_enc <- max(abs(enc_direct - enc_mizer) / enc_mizer)
fl <- getFeedingLevel(params)
Q <- matrix(0, no_sp, no_w_full)
Q[, idx_sp] <- (1 - fl) * params@search_vol * sweep(n, 2, dw(params), "*")
pr_direct <- t(sapply(seq_len(no_sp), function(i) circular(phi_p, Q[i, ], i, seq_len(no_w_full))))
pr_direct[pr_direct < 1e-18] <- 0
pr_direct <- pr_direct * params@ft_mask
pr_mizer <- getPredRate(params)
rel_err_pr <- max(abs(pr_direct - pr_mizer)) / max(pr_mizer)
cat(sprintf("Circular sums vs mizer FFT: encounter %.2e, pred rate %.2e\n",
            rel_err_enc, rel_err_pr))
stopifnot(rel_err_enc < 1e-10, rel_err_pr < 1e-10)

# Unfished states, derived exactly as in the paper
states <- get_fixed_points(params)
state_json <- lapply(states, function(p) {
  list(n = unname(initialN(p)), n_pp = unname(as.numeric(initialNResource(p))),
       biomass = unname(as.numeric(getBiomass(p))))
})

# Fishing mortality per unit effort for each gear, as species x size
sel <- params@selectivity       # gear x species x size
F_per_effort <- matrix(0, no_sp, no_w)
gear_of_species <- integer(no_sp)
for (g in seq_len(dim(sel)[1])) {
  for (i in seq_len(no_sp)) {
    s <- sel[g, i, ] * params@catchability[g, i]
    if (any(s > 0)) gear_of_species[i] <- g
    F_per_effort[i, ] <- F_per_effort[i, ] + s
  }
}
stopifnot(all(gear_of_species == seq_len(no_sp)))  # one gear per species

model <- list(
  species = sp$species,
  colours = unname(getColours(params)[sp$species]),
  w = w(params), dw = dw(params),
  w_full = w_full(params), dw_full = dw_full(params),
  offset0 = offset0,
  w_min_idx = unname(params@w_min_idx) - 1,       # 0-based
  w_mat = sp$w_mat, w_max = sp$w_max, w_min = sp$w_min,
  w_max_idx = unname(vapply(sp$w_max, function(wm) sum(w(params) <= wm), 1L)),
  alpha = sp$alpha, erepro = sp$erepro, R_max = sp$R_max,
  interaction = unname(interaction_matrix(params)),
  interaction_resource = sp$interaction_resource,
  search_vol = unname(params@search_vol),
  intake_max = unname(params@intake_max),
  metab = unname(params@metab),
  psi = unname(params@psi),
  mu_b = unname(params@mu_b),
  F_per_effort = unname(F_per_effort),
  knife_edge_size = gear_params(params)$knife_edge_size,
  phi_e = unname(phi_e),
  phi_p = unname(phi_p),
  ft_mask = unname(params@ft_mask * 1),
  rr_pp = as.numeric(params@rr_pp),
  cc_pp = as.numeric(params@cc_pp),
  states = state_json
)

writeLines(
  c("// Generated by interactive/export_model.R from cultivation_depensation_params.rds.",
    "// Do not edit by hand.",
    paste0("(function (root) { var MODEL = ",
           toJSON(model, digits = I(15), auto_unbox = TRUE, matrix = "rowmajor"),
           "; if (typeof module !== 'undefined') module.exports = MODEL;",
           " else root.MODEL = MODEL; })(this);")),
  "interactive/model_data.js"
)

# ---------------------------------------------------------------------------
# Reference trajectories for validating the JavaScript engine
# ---------------------------------------------------------------------------
biomass_ts <- function(sim) {
  b <- getBiomass(sim)
  list(t = as.numeric(rownames(b)), forage = unname(b[, "Forage"]),
       predator = unname(b[, "Predator"]))
}
pA <- states$State_A
pB <- states$State_B
pC <- states$State_C

collapse <- project(pA, effort = c(Forage = 0, Predator = 0.25), t_max = 150,
                    dt = 0.2, t_save = 1, progress_bar = FALSE)
balanced <- project(pA, effort = c(Forage = 0.35, Predator = 0.25), t_max = 150,
                    dt = 0.2, t_save = 1, progress_bar = FALSE)
inv_B <- pB
initialN(inv_B)["Forage", ] <- initialN(pA)["Forage", ] * 1e-3
invade_B <- project(inv_B, t_max = 150, dt = 0.2, t_save = 1, progress_bar = FALSE)
res_C <- pC
initialN(res_C)["Predator", ] <- initialN(pA)["Predator", ] * 5e-3
rescue_C <- project(res_C, effort = c(Forage = 0.30, Predator = 0), t_max = 150,
                    dt = 0.2, t_save = 1, progress_bar = FALSE)

r0 <- lapply(states, function(p) as.list(get_R0(p)))
r0_C_fished <- as.list(get_R0(finalParams(
  project(pC, effort = c(Forage = 0.30, Predator = 0), t_max = 300, dt = 0.2,
          progress_bar = FALSE)),
  effort = c(Forage = 0.30, Predator = 0)))

reference <- list(
  mizer_version = as.character(packageVersion("mizer")),
  dt = 0.2,
  collapse = c(biomass_ts(collapse), effort = list(c(0, 0.25)), start = "State_A"),
  balanced = c(biomass_ts(balanced), effort = list(c(0.35, 0.25)), start = "State_A"),
  invade_B = c(biomass_ts(invade_B), effort = list(c(0, 0)), start = "State_B",
               seed = list(list(species = 0, factor = 1e-3))),
  rescue_C = c(biomass_ts(rescue_C), effort = list(c(0.30, 0)), start = "State_C",
               seed = list(list(species = 1, factor = 5e-3))),
  R0 = r0,
  R0_C_forage_0.30 = r0_C_fished,
  rates_A = list(
    encounter = unname(getEncounter(pA)), e_growth = unname(getEGrowth(pA)),
    mort = unname(getMort(pA)), rdd = unname(getRDD(pA)),
    resource_mort = unname(as.numeric(getResourceMort(pA)))
  )
)
write_json(reference, "interactive/validation/reference.json", digits = I(15),
           auto_unbox = TRUE, matrix = "rowmajor")
cat("Wrote interactive/model_data.js and interactive/validation/reference.json\n")
print(sapply(states, function(p) round(getBiomass(p), 2)))
