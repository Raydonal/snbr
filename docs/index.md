# snbr

Bias reduction in linear and partially linear regression models with
skew-normal errors.

The maximum likelihood estimate of the shape parameter of the
skew-normal distribution is infinite with positive probability in finite
samples, and when it is finite its bias can be large. `snbr` provides

  - maximum likelihood estimation with an exact test of whether the
    estimate of the shape is infinite (comparison with the limiting
    half-normal regression, computed by constrained least squares);
  - the second-order (Cox and Snell) bias of all parameters, computed
    from exact cumulants of the log-likelihood derivatives;
  - the corrective estimator and the mean (Firth) and median (Kenne
    Pagui, Salvan and Sartori) bias-reducing adjusted score estimators,
    for the location coefficients, scale and shape jointly;
  - a partially linear model with a P-spline term, with the second-order
    bias of the penalized estimator at a fixed smoothing parameter.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("raydonal/snbr")
```

## Example

``` r
library(snbr)
arcticlake$d <- arcticlake$depth / 100
snreg(logit_clay ~ d + I(d^2), data = arcticlake, method = "ML")   # -Inf
snreg(logit_clay ~ d + I(d^2), data = arcticlake, method = "mean") # finite
snreg(logit_clay ~ 1, data = arcticlake, smooth = "d", method = "mean")
```

## Reference

Lopes, E. A. C., Vasconcellos, K. L. P. and Ospina, R. Bias reduction
and infinite shape estimates in linear and partially linear skew-normal
regression. Manuscript.
