# Loaded by paper.Rmd: every number quoted in the paper comes from the results
# file written by cultivation_depensation.R.
res <- readRDS("results/cultivation_results.rds")

# Fixed number of decimals
fmt <- function(x, digits = 2) formatC(unname(x), format = "f", digits = digits)
# Significant figures
sig <- function(x, digits = 2) format(signif(unname(x), digits), big.mark = ",")

state <- function(nm) res$states[res$states$State == nm, ]
bound <- function(Ff) res$bounds[abs(res$bounds$F_Forage - Ff) < 1e-9, ]
r0C <- function(Ff) res$r0_curve$R0_Pred_C[abs(res$r0_curve$F_Forage - Ff) < 1e-9]
BC <- function(Ff) res$r0_curve$Forage_Biomass_C[abs(res$r0_curve$F_Forage - Ff) < 1e-9]
sp <- function(col, species) res$species_params[res$species_params$species == species, col]
# The hysteresis closes between the last grid point with a bistable interval and
# the first without one
cusp_range <- sprintf("%s--%s", fmt(res$cusp_F_Forage - 0.025, 3), fmt(res$cusp_F_Forage, 3))
hyst_width <- function(Ff) bound(Ff)$F_up - bound(Ff)$F_low
# Size diffusion of the reference model (the default of build_cultivation_model())
if (is.null(res$c_diff)) res$c_diff <- 0.1
