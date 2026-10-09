# Partially linear skew-normal regression with bias reduction

Fits \\(y\_i = x\_i^\\top\\beta + f(t\_i) + \\varepsilon\_i\\),
\\(\\varepsilon\_i \\sim SN(0, \\sigma, \\alpha)\\), with \\(f\\) a
P-spline. For a fixed smoothing parameter the estimating equation is the
penalized score, optionally corrected by the mean bias-reducing
adjustment computed with the penalized information \\(H = K + \\lambda
S\\). The smoothing parameter is chosen on a grid by AIC with effective
degrees of freedom \\(\\mathrm{tr}(H^{-1}K)\\).

## Usage

``` r
snplm(
  X,
  t,
  y,
  method = c("BR", "ML"),
  k = 12L,
  lambda = NULL,
  lambda_grid = 10^seq(-3, 4, by = 0.5)
)
```

## Arguments

  - X:
    
    linear part of the model matrix (including the intercept).

  - t:
    
    covariate of the smooth term.

  - y:
    
    response.

  - method:
    
    \`"ML"\` (penalized likelihood) or \`"BR"\` (penalized likelihood
    with the mean bias-reducing adjustment).

  - k:
    
    basis dimension.

  - lambda:
    
    optional fixed smoothing parameter; if \`NULL\` it is selected from
    \`lambda\_grid\`.

  - lambda\_grid:
    
    grid of candidate values.

## Value

A list with the coefficients, the smoothing parameter, the effective
degrees of freedom, the fitted smooth, AIC and the basis.
