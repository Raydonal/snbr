# snbr 0.2.1

* Bug fix: `sn_bias()` assigned wrong/NA parameter names (e.g. `sigma`,
  `alpha`, `NA`, `NA`) whenever the model matrix had no column names,
  because it indexed a too-short name vector with `seq_along()`. `sn_info()`
  and the matrix interface of `sn_mle()`, `sn_br()` and `snplm()` had the
  milder version of the same bug (silently unnamed output instead of wrong
  names). Also fixed: `cbind(1, d, d^2)` sets `colnames()` to `c("", "d",
  "")`, not `NULL` (R auto-names columns built from a bare symbol and
  leaves literal/expression columns blank), which defeated the previous
  `is.null(colnames(X))` fallback checks. All five sites now go through a
  single internal helper (`.coef_names()`, `R/utils.R`) that fills in
  `beta1`, `beta2`, ... at every missing or blank position. No numeric
  result changes: verified against `analises/R/01_verifica_formula_artigo.R`
  and `R/13_aplicacao_lago.R`, which reproduce the same values as before,
  now with correct labels.
* `print()`/`summary()` on an `"snbr"` object now explain, when the maximum
  likelihood estimate of the shape is infinite, that this is an exact
  result rather than a numerical failure, and suggest refitting with
  `method = "mean"` or `"median"`.

# snbr 0.2.0

* Version that accompanies the manuscript submitted to *Information and
  Inference*; every result in the manuscript and its supplement was produced
  with this version by the scripts in `analises/R`.
* `sn_br()` reports why an adjusted score has no root (`failure`): a sign
  change only across alpha = 0, no sign change, or nonconvergence.
* `inst/CITATION` and documentation updated to the submitted manuscript.

# snbr 0.1.0

* First version: maximum likelihood with exact detection of infinite shape
  estimates, Cox and Snell bias, corrective estimator, mean and median
  bias-reducing adjusted scores, parametric bootstrap correction, and a
  partially linear model with a P-spline term.
