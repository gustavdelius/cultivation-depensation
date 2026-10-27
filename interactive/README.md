# The Nursery Trap

A playable version of the cultivation–depensation model in `paper.Rmd`, aimed at
fishery policy makers. Open `nursery_trap.html` in any browser; it is a single
self-contained file that works offline and can be emailed.

The page runs the paper's mizer model live in the browser. Nothing is scripted:
the stocks, catches, survival rates and comeback number all come from the model
as the player moves the fishing levers.

## Files

| File | Role |
|---|---|
| `app.html` | Page source (markup, styles, game logic) |
| `engine.js` | JavaScript port of mizer's default Euler/upwind dynamics |
| `model_data.js` | Rate arrays and the unfished states 0, A, B, C, exported from mizer |
| `regimes_data.js` | Collapse and invasion boundaries and equilibrium yields, computed with the engine |
| `nursery_trap.html` | Built, self-contained page (generated) |
| `export_model.R` | Writes `model_data.js` and `validation/reference.json` from `cultivation_depensation_params.rds` |
| `validate.js` | Checks the engine against mizer trajectories, rates and R0 values |
| `compute_regimes.js` | Writes `regimes_data.js` |
| `build.js` | Inlines the scripts into `app.html` to make `nursery_trap.html` |

## Rebuilding

Run from the project root:

```bash
Rscript interactive/export_model.R
node interactive/validate.js
node interactive/compute_regimes.js
node interactive/build.js
```

`validate.js` must report "All checks passed": the engine reproduces mizer's
rates, R0 values and 150-year trajectories to a relative error below 1e-13.

## Faithfulness to mizer

mizer evaluates encounter and predation rates by FFT, which makes them circular
convolutions over the `w_full` grid. The engine evaluates the same circular sums
directly, wrap-around included, so that it matches mizer exactly. In this model
the wrap-around adds about 9% to the encounter rate of 1 mg forage larvae at
State A, because the resource grid starts only three decades below the egg size.
