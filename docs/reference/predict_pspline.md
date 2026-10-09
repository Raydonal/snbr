# Evaluate a fitted P-spline component at new points

Evaluate a fitted P-spline component at new points

## Usage

``` r
predict_pspline(basis, gamma, tnew)
```

## Arguments

  - basis:
    
    output of \[pspline\_basis()\].

  - gamma:
    
    coefficients of the constrained basis.

  - tnew:
    
    new values of the covariate.

## Value

Vector of fitted values of \\(f\\) at \`tnew\`.
