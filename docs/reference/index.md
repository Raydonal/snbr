# Package index

## Model fitting

Maximum likelihood with exact detection of infinite shape estimates, and
the corrective and preventive bias-reduced estimators.

- [`snreg()`](https://raydonal.github.io/snbr/reference/snreg.md) :
  Skew-normal regression with bias reduction
- [`sn_mle()`](https://raydonal.github.io/snbr/reference/sn_mle.md) :
  Maximum likelihood fit of the skew-normal regression model
- [`sn_br()`](https://raydonal.github.io/snbr/reference/sn_br.md) :
  Bias-reducing adjusted score estimators
- [`sn_correct()`](https://raydonal.github.io/snbr/reference/sn_correct.md)
  : Corrective (Cox and Snell) bias-corrected estimator
- [`sn_boot()`](https://raydonal.github.io/snbr/reference/sn_boot.md) :
  Parametric bootstrap bias correction
- [`snplm()`](https://raydonal.github.io/snbr/reference/snplm.md) :
  Partially linear skew-normal regression with bias reduction

## Second-order bias and information

- [`sn_bias()`](https://raydonal.github.io/snbr/reference/sn_bias.md) :
  Second-order bias of the (penalized) maximum likelihood estimator
- [`sn_info()`](https://raydonal.github.io/snbr/reference/sn_info.md) :
  Expected information of the skew-normal regression model
- [`sn_cumulants()`](https://raydonal.github.io/snbr/reference/sn_cumulants.md)
  : Exact unit cumulants of the skew-normal log-likelihood derivatives
- [`sn_cumulants_fast()`](https://raydonal.github.io/snbr/reference/sn_cumulants_fast.md)
  : Unit cumulants by interpolation

## Likelihood

- [`sn_loglik()`](https://raydonal.github.io/snbr/reference/sn_loglik.md)
  [`sn_score()`](https://raydonal.github.io/snbr/reference/sn_loglik.md)
  [`sn_hessian()`](https://raydonal.github.io/snbr/reference/sn_loglik.md)
  : Log-likelihood and score of the skew-normal regression model
- [`rsnorm()`](https://raydonal.github.io/snbr/reference/rsnorm.md) :
  Random generation from the skew-normal distribution
- [`psnorm()`](https://raydonal.github.io/snbr/reference/psnorm.md) :
  Skew-normal distribution function

## Smooth terms

- [`pspline_basis()`](https://raydonal.github.io/snbr/reference/pspline_basis.md)
  : P-spline basis with a sum-to-zero constraint
- [`predict_pspline()`](https://raydonal.github.io/snbr/reference/predict_pspline.md)
  : Evaluate a fitted P-spline component at new points

## Data

- [`arcticlake`](https://raydonal.github.io/snbr/reference/arcticlake.md)
  : Arctic lake sediments
