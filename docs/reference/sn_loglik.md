# Log-likelihood and score of the skew-normal regression model

Log-likelihood and score of the skew-normal regression model

## Usage

``` r
sn_loglik(theta, X, y)

sn_score(theta, X, y)

sn_hessian(theta, X, y)
```

## Arguments

  - theta:
    
    \`c(beta, sigma, alpha)\`.

  - X:
    
    model matrix (n x p).

  - y:
    
    response vector.

## Value

\`sn\_loglik()\` returns a scalar, \`sn\_score()\` a vector of length
\\(p+2\\), \`sn\_hessian()\` the observed Hessian matrix.
