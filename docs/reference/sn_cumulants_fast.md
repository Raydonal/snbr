# Unit cumulants by interpolation

Same output as \[sn_cumulants()\], obtained by spline interpolation of a
precomputed grid. The interpolation error is below \\10^{-8}\\ in
relative terms over \\\|\alpha\| \le 2000\\ (checked in the package
tests); outside that range the exact quadrature is used.

## Usage

``` r
sn_cumulants_fast(alpha, sigma = 1)
```

## Arguments

- alpha:

  shape parameter (scalar).

- sigma:

  scale parameter (scalar, default 1).

## Value

See \[sn_cumulants()\].
