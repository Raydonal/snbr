# Maximum likelihood fit of the skew-normal regression model

Maximizes the log-likelihood from several starting points and then
decides whether the shape estimate is finite. The decision is exact: as
\\\alpha \to \pm\infty\\ the log-likelihood converges to that of a
half-normal regression whose residuals all have the sign of \\\alpha\\,
and its supremum is computed by constrained least squares. When this
limit is not smaller than the best finite value, the maximum likelihood
estimate of \\\alpha\\ is reported as \`Inf\` or \`-Inf\`.

## Usage

``` r
sn_mle(X, y, start = NULL, tol = 1e-08)
```

## Arguments

- X:

  model matrix (n x p).

- y:

  response.

- start:

  optional starting value \`c(beta, sigma, alpha)\`.

- tol:

  tolerance (on the log-likelihood scale) used in the comparison with
  the limit.

## Value

A list with \`theta\`, \`loglik\`, \`finite\` (logical) and
\`loglik_limit\` (the two limits, for \`alpha -\> -Inf, +Inf\`).
