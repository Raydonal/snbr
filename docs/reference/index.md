# Package index

## Model fitting

Maximum likelihood with exact detection of infinite shape estimates, and
the corrective and preventive bias-reduced estimators.

<!-- end list -->

  - `snreg()` : Skew-normal regression with bias reduction
  - `sn_mle()` : Maximum likelihood fit of the skew-normal regression
    model
  - `sn_br()` : Bias-reducing adjusted score estimators
  - `sn_correct()` : Corrective (Cox and Snell) bias-corrected estimator
  - `sn_boot()` : Parametric bootstrap bias correction
  - `snplm()` : Partially linear skew-normal regression with bias
    reduction

## Second-order bias and information

<!-- end list -->

  - `sn_bias()` : Second-order bias of the (penalized) maximum
    likelihood estimator
  - `sn_info()` : Expected information of the skew-normal regression
    model
  - `sn_cumulants()` : Exact unit cumulants of the skew-normal
    log-likelihood derivatives
  - `sn_cumulants_fast()` : Unit cumulants by interpolation

## Likelihood

<!-- end list -->

  - `sn_loglik()` `sn_score()` `sn_hessian()` : Log-likelihood and score
    of the skew-normal regression model
  - `rsnorm()` : Random generation from the skew-normal distribution
  - `psnorm()` : Skew-normal distribution function

## Smooth terms

<!-- end list -->

  - `pspline_basis()` : P-spline basis with a sum-to-zero constraint
  - `predict_pspline()` : Evaluate a fitted P-spline component at new
    points

## Data

<!-- end list -->

  - `arcticlake` : Arctic lake sediments
