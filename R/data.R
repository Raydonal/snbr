#' Arctic lake sediments
#'
#' Sand, silt and clay proportions of 39 sediment samples collected at
#' different water depths in an Arctic lake (Aitchison, 1986). The clay
#' proportion on the logit scale is the response of the illustrative
#' example in the package vignette.
#'
#' @format A data frame with 39 rows and 5 variables:
#' \describe{
#'   \item{depth}{water depth (m).}
#'   \item{sand, silt, clay}{proportions of the three grain sizes.}
#'   \item{logit_clay}{\code{log(clay / (1 - clay))}.}
#' }
#' @source Aitchison, J. (1986). *The Statistical Analysis of Compositional
#'   Data*. Chapman and Hall, London. The values coincide with the data set
#'   `ArcticLake` of the package compositions.
"arcticlake"
