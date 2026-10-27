# Loaded by paper.Rmd: every number quoted in the paper comes from the results
# file written by cultivation_depensation.R.
res <- readRDS("results/cultivation_results.rds")

# Fixed number of decimals
fmt <- function(x, digits = 2) formatC(unname(x), format = "f", digits = digits)
# Significant figures
sig <- function(x, digits = 2) format(signif(unname(x), digits), big.mark = ",")
# Scientific notation for LaTeX math
sci <- function(x, digits = 2) {
  e <- floor(log10(abs(unname(x))))
  sprintf("%s \\times 10^{%d}", formatC(unname(x) / 10^e, format = "f", digits = digits - 1), e)
}
# Power of ten for LaTeX math
pow10 <- function(x) sprintf("10^{%d}", round(log10(unname(x))))

state <- function(nm) res$states[res$states$State == nm, ]
bound <- function(Ff) res$bounds[abs(res$bounds$F_Forage - Ff) < 1e-9, ]
r0C <- function(Ff) res$r0_curve$R0_Pred_C[abs(res$r0_curve$F_Forage - Ff) < 1e-9]
BC <- function(Ff) res$r0_curve$Forage_Biomass_C[abs(res$r0_curve$F_Forage - Ff) < 1e-9]
sp <- function(col, species) res$species_params[res$species_params$species == species, col]
# The hysteresis closes between the last grid point with a bistable interval and
# the first without one
cusp_range <- sprintf("%s\\text{--}%s", fmt(res$cusp_F_Forage - 0.025, 3), fmt(res$cusp_F_Forage, 3))
hyst_width <- function(Ff) bound(Ff)$F_up - bound(Ff)$F_low
# F_Forage above which the invasion boundary of State C lies above the fold of
# the abundant state, so that its collapse leads to the low-predator state
inv_above_up <- with(res$bounds[!is.na(res$bounds$F_inv), ], {
  d <- F_inv - F_up; i <- which(d > 0)[1]
  if (is.na(i)) NA else if (i == 1) F_Forage[1] else approx(d[(i - 1):i], F_Forage[(i - 1):i], xout = 0)$y
})
# Size diffusion of the reference model (the default of build_cultivation_model())
if (is.null(res$c_diff)) res$c_diff <- 0.1
