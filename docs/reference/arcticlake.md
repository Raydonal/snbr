# Arctic lake sediments

Sand, silt and clay proportions of 39 sediment samples collected at
different water depths in an Arctic lake (Aitchison, 1986). The clay
proportion on the logit scale is the response of the illustrative
example in the package vignette.

## Usage

``` r
arcticlake
```

## Format

A data frame with 39 rows and 5 variables:

- depth:

  water depth (m).

- sand, silt, clay:

  proportions of the three grain sizes.

- logit_clay:

  `log(clay / (1 - clay))`.

## Source

Aitchison, J. (1986). \*The Statistical Analysis of Compositional
Data\*. Chapman and Hall, London. The values coincide with the data set
\`ArcticLake\` of the package compositions.
