# P-spline basis with a sum-to-zero constraint

P-spline basis with a sum-to-zero constraint

## Usage

``` r
pspline_basis(t, k = 12L, degree = 3L, diff = 2L)
```

## Arguments

- t:

  covariate entering nonparametrically.

- k:

  number of B-spline basis functions before the constraint.

- degree:

  spline degree.

- diff:

  order of the difference penalty.

## Value

A list with the constrained basis \`Z\` (n x (k-1)), the penalty \`S\`
((k-1) x (k-1)) and the objects needed by \[predict_pspline()\].

## References

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties. \*Statistical Science\*, 11, 89–121.
