# Corrective (Cox and Snell) bias-corrected estimator

Returns \\(\\tilde\\theta = \\hat\\theta - b(\\hat\\theta)\\), with
\\(b\\) the second-order bias of \[sn\_bias()\]. Not defined when the
maximum likelihood estimate of the shape is infinite.

## Usage

``` r
sn_correct(theta, X)
```

## Arguments

  - theta:
    
    maximum likelihood estimate.

  - X:
    
    model matrix.

## Value

Numeric vector (all \`NA\` if \`theta\` is not finite).
