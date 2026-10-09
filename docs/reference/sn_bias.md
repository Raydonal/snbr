# Second-order bias of the (penalized) maximum likelihood estimator

Evaluates the \\(O(n^{-1})\\) bias of the maximum likelihood estimator
of \\(\\theta = (\\beta, \\sigma, \\alpha)\\) in the skew-normal
regression model (Cox and Snell, 1968). When a quadratic penalty
\\(\\lambda \\theta^\\top S \\theta / 2\\) is present, the bias is that
of the penalized estimator around its own target, obtained from the same
expansion with the penalized information \\(H = K + \\lambda S\\).

## Usage

``` r
sn_bias(theta, X, S = NULL, lambda = 0, cu = NULL)
```

## Arguments

  - theta:
    
    parameter vector \`c(beta, sigma, alpha)\`.

  - X:
    
    model matrix of the location (n x p).

  - S:
    
    optional penalty matrix of dimension \\((p+2) \\times (p+2)\\) (zero
    rows and columns for unpenalized parameters).

  - lambda:
    
    smoothing parameter (non-negative scalar).

  - cu:
    
    optional output of \[sn\_cumulants()\] at \`theta\`'s shape and
    scale, to avoid recomputation.

## Value

Numeric vector of length \\(p+2\\) with the estimated biases.

## References

Cox, D. R. and Snell, E. J. (1968). A general definition of residuals.
\*Journal of the Royal Statistical Society B\*, 30, 248–275.
