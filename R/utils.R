## Fills in a generic name at every NULL/NA/"" position of `nm` (length p),
## leaving any name the caller actually supplied untouched. Needed because
## cbind(1, d, d^2) sets colnames to c("", "d", "") rather than NULL: it
## auto-names columns built from a bare symbol and leaves literal/expression
## columns blank, so a plain is.null(colnames(X)) check misses this case.
.fill_blanks <- function(nm, p, prefix) {
  if (is.null(nm)) nm <- rep.int("", p)
  blank <- is.na(nm) | !nzchar(nm)
  if (any(blank)) nm[blank] <- paste0(prefix, seq_len(p))[blank]
  nm
}

## Coefficient names for (beta, [gamma,] sigma, alpha). Falls back to
## generic names, column by column, when a model matrix has no (or no
## usable) column names, so that callers of the matrix interface
## (sn_mle(), sn_bias(), sn_br(), snplm(), ...) on a bare matrix never get
## silently misaligned, blank or NA names.
.coef_names <- function(X, extra = NULL, extra_prefix = "g") {
  nm <- .fill_blanks(colnames(X), ncol(X), "beta")
  if (!is.null(extra)) {
    nm <- c(nm, .fill_blanks(colnames(extra), ncol(extra), extra_prefix))
  }
  c(nm, "sigma", "alpha")
}
