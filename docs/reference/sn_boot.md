# Parametric bootstrap bias correction

Generates \`B\` samples from the model fitted by maximum likelihood,
refits, and returns \\(2\\hat\\theta - \\bar\\theta^\*\\). Bootstrap
replicates whose shape estimate is infinite are counted and, following
the policy chosen by \`infinite\`, either discarded (the policy of the
earlier versions of this work) or make the correction undefined.

## Usage

``` r
sn_boot(theta, X, B = 500L, infinite = c("discard", "na"))
```

## Arguments

  - theta:
    
    maximum likelihood estimate (finite).

  - X:
    
    model matrix.

  - B:
    
    number of bootstrap samples.

  - infinite:
    
    \`"discard"\` or \`"na"\`.

## Value

A list with the corrected estimate \`theta\`, the number of infinite
replicates \`n\_inf\` and the bootstrap mean \`mean\_star\`.
