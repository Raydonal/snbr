# Bias-reducing adjusted score estimators

Solves the adjusted score equations of Firth (1993) (\`type = "mean"\`)
or of Kenne Pagui, Salvan and Sartori (2017) (\`type = "median"\`) for
all parameters of the skew-normal regression model, by the quasi
Fisher-scoring iteration of Kosmidis and Firth (2010). A quadratic
penalty \\(\\lambda \\theta^\\top S \\theta/2\\) can be added (mean
adjustment only), as used for the semiparametric model.

## Usage

``` r
sn_br(
  X,
  y,
  type = c("mean", "median"),
  start = NULL,
  S = NULL,
  lambda = 0,
  maxit = 200L,
  tol = 1e-08
)
```

## Arguments

  - X:
    
    model matrix.

  - y:
    
    response.

  - type:
    
    \`"mean"\` or \`"median"\`.

  - start:
    
    optional starting value.

  - S, lambda:
    
    optional penalty.

  - maxit, tol:
    
    iteration controls.

## Value

A list with \`theta\`, \`converged\`, \`iterations\` and the value of
the adjusted score at the solution.

## References

Firth, D. (1993). Bias reduction of maximum likelihood estimates.
\*Biometrika\*, 80, 27–38.

Kosmidis, I. and Firth, D. (2010). A generic algorithm for reducing bias
in parametric estimation. \*Electronic Journal of Statistics\*, 4,
1097–1112.

Kenne Pagui, E. C., Salvan, A. and Sartori, N. (2017). Median bias
reduction of maximum likelihood estimates. \*Biometrika\*, 104, 923–938.
