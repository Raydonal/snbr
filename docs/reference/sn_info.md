# Expected information of the skew-normal regression model

Expected information of the skew-normal regression model

## Usage

``` r
sn_info(theta, X, cu = NULL)
```

## Arguments

  - theta:
    
    parameter vector \`c(beta, sigma, alpha)\`.

  - X:
    
    model matrix of the location (n x p).

  - cu:
    
    optional output of \[sn\_cumulants()\] at \`theta\`'s shape and
    scale, to avoid recomputation.

## Value

The \\((p+2) \\times (p+2)\\) Fisher information matrix.
