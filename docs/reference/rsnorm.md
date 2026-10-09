# Random generation from the skew-normal distribution

Uses the stochastic representation \\(Z = \\delta |U| +
(1-\\delta^2)^{1/2} V\\), \\(\\delta = \\alpha/(1+\\alpha^2)^{1/2}\\),
with \\(U, V\\) independent standard normal (Henze, 1986).

## Usage

``` r
rsnorm(n, xi = 0, omega = 1, alpha = 0)
```

## Arguments

  - n:
    
    number of observations.

  - xi, omega, alpha:
    
    location, scale and shape (recycled).

## Value

Numeric vector of length \`n\`.

## Examples

``` r
z <- rsnorm(1000, alpha = 5)
mean(z)  # close to sqrt(2/pi) * 5 / sqrt(26)
#> [1] 0.7845761
```
