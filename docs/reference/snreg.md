# Skew-normal regression with bias reduction

Fits the linear model \\(y\_i = x\_i^\\top\\beta + \\varepsilon\_i\\) or
the partially linear model \\(y\_i = x\_i^\\top\\beta + f(t\_i) +
\\varepsilon\_i\\), with \\(\\varepsilon\_i \\sim SN(0, \\sigma,
\\alpha)\\) independent, by one of four methods:

  - \`"ML"\`:
    
    maximum likelihood; the shape estimate is reported as \`Inf\` or
    \`-Inf\` when the likelihood has no finite maximizer (see
    \[sn\_mle()\]).

  - \`"CS"\`:
    
    the corrective estimator \\(\\hat\\theta - b(\\hat\\theta)\\) with
    the Cox and Snell bias \\(b\\) (\[sn\_bias()\]).

  - \`"mean"\`:
    
    the mean bias-reducing adjusted score estimator (Firth, 1993),
    always finite in our experience.

  - \`"median"\`:
    
    the median bias-reducing adjusted score estimator (Kenne Pagui,
    Salvan and Sartori, 2017); linear models only.

With a smooth term only \`"ML"\` (penalized likelihood) and \`"mean"\`
are available; see \[snplm()\].

## Usage

``` r
snreg(
  formula,
  data,
  smooth = NULL,
  method = c("mean", "median", "ML", "CS"),
  k = 10L,
  lambda = NULL
)
```

## Arguments

  - formula:
    
    model formula for the location (linear part).

  - data:
    
    data frame.

  - smooth:
    
    optional name (character) of a covariate entering through a P-spline
    term.

  - method:
    
    one of \`"mean"\`, \`"median"\`, \`"ML"\`, \`"CS"\`.

  - k:
    
    basis dimension of the smooth term.

  - lambda:
    
    optional fixed smoothing parameter.

## Value

An object of class \`"snbr"\`.

## Examples

``` r
fit <- snreg(logit_clay ~ depth, data = arcticlake, method = "mean")
summary(fit)
#> snreg(formula = logit_clay ~ depth, data = arcticlake, method = "mean")
#> Method: mean 
#> 
#>               Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) -1.3796040  0.3311546 -4.1660 3.099e-05 ***
#> depth        0.0282985  0.0048912  5.7856 7.223e-09 ***
#> sigma        1.5473465  0.2502512                      
#> alpha       -3.4151515  1.9385987 -1.7617   0.07813 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -49.5592  (effective df 4 )
```
