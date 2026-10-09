## Builds the interpolation grid used by sn_cumulants_fast().
## Exact unit cumulants (sigma = 1) are computed by adaptive quadrature at
## alpha = sinh(t), t on a regular grid in [0, asinh(2000)], and stored in
## R/sysdata.rda. Negative alpha is handled by the reflection
## k(-alpha) = (-1)^{#eta + #alpha indices} k(alpha).
## Run from the package root:  Rscript data-raw/cumulant_grid.R
source("R/derivs.R")
source("R/cumulants.R")
tgrid <- seq(0, asinh(2000), length.out = 1601)
agrid <- sinh(tgrid)
vals <- parallel::mclapply(agrid, function(a) {
  cu <- sn_cumulants(a)
  c(as.vector(cu$k2), as.vector(cu$k3), as.vector(cu$k21), as.vector(cu$k111))
}, mc.cores = 4)
vals <- do.call(rbind, vals)
.cumgrid <- list(t = tgrid, vals = vals, built = Sys.time(),
                 note = "columns: k2 (9), k3 (27), k21 (27), k111 (27), column-major")
save(.cumgrid, file = "R/sysdata.rda", compress = "xz")
cat("grid saved:", nrow(vals), "x", ncol(vals), "\n")
