# Exact unit cumulants of the skew-normal log-likelihood derivatives

Computes, by adaptive quadrature, the expected values of products of
derivatives of the skew-normal log-density with respect to the reduced
parameter \\(\eta, \sigma, \alpha)\\, where \\\eta\\ is the location.
These are the building blocks of the Fisher information, of the Cox and
Snell second-order bias and of the bias-reducing adjusted scores.

## Usage

``` r
sn_cumulants(alpha, sigma = 1, rel.tol = 1e-11)
```

## Arguments

- alpha:

  shape parameter (scalar).

- sigma:

  scale parameter (scalar, default 1).

- rel.tol:

  relative tolerance passed to \[stats::integrate()\].

## Value

A list with components \`k2\` (3 x 3), \`k3\`, \`k21\` and \`k111\` (3 x
3 x 3 arrays), indexed by \`c("eta", "sigma", "alpha")\`. For \`k21\`
the first two indices refer to the second derivative and the third to
the first derivative.

## Examples

``` r
cu <- sn_cumulants(3)
-cu$k2   # unit Fisher information at sigma = 1
#>             eta       sigma       alpha
#> eta   3.0281308  1.03041365  0.16115529
#> sigma 1.0304137  2.24420447 -0.08140149
#> alpha 0.1611553 -0.08140149  0.02713383
```
