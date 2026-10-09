#' Skew-normal regression with bias reduction
#'
#' Fits the linear model \eqn{y_i = x_i^\top\beta + \varepsilon_i} or the
#' partially linear model \eqn{y_i = x_i^\top\beta + f(t_i) +
#' \varepsilon_i}, with \eqn{\varepsilon_i \sim SN(0, \sigma, \alpha)}
#' independent, by one of four methods:
#' \describe{
#'   \item{`"ML"`}{maximum likelihood; the shape estimate is reported as
#'     `Inf` or `-Inf` when the likelihood has no finite maximizer
#'     (see [sn_mle()]).}
#'   \item{`"CS"`}{the corrective estimator \eqn{\hat\theta -
#'     b(\hat\theta)} with the Cox and Snell bias \eqn{b} ([sn_bias()]).}
#'   \item{`"mean"`}{the mean bias-reducing adjusted score estimator
#'     (Firth, 1993), always finite in our experience.}
#'   \item{`"median"`}{the median bias-reducing adjusted score estimator
#'     (Kenne Pagui, Salvan and Sartori, 2017); linear models only.}
#' }
#' With a smooth term only `"ML"` (penalized likelihood) and `"mean"` are
#' available; see [snplm()].
#'
#' @param formula model formula for the location (linear part).
#' @param data data frame.
#' @param smooth optional name (character) of a covariate entering
#'   through a P-spline term.
#' @param method one of `"mean"`, `"median"`, `"ML"`, `"CS"`.
#' @param k basis dimension of the smooth term.
#' @param lambda optional fixed smoothing parameter.
#' @return An object of class `"snbr"`.
#' @examples
#' fit <- snreg(logit_clay ~ depth, data = arcticlake, method = "mean")
#' summary(fit)
#' @export
snreg <- function(formula, data, smooth = NULL, method = c("mean", "median", "ML", "CS"),
                  k = 10L, lambda = NULL) {
  method <- match.arg(method)
  cl <- match.call()
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, mf)
  smooth_info <- NULL
  finite <- NA
  if (!is.null(smooth)) {
    if (!method %in% c("mean", "ML")) stop("with a smooth term use method = \"mean\" or \"ML\"")
    t <- data[[smooth]]
    fit <- snplm(X, t, y, method = if (method == "mean") "BR" else "ML", k = k, lambda = lambda)
    theta <- fit$theta
    Xf <- fit$X
    S <- fit$S
    lam <- fit$lambda
    edf <- fit$edf
    smooth_info <- list(var = smooth, basis = fit$basis, gamma = fit$gamma, t = t)
  } else {
    Xf <- X
    S <- NULL
    lam <- 0
    if (method %in% c("ML", "CS")) {
      f <- sn_mle(X, y)
      finite <- f$finite
      theta <- if (method == "ML") f$theta else sn_correct(f$theta, X)
    } else {
      f <- sn_br(X, y, method)
      if (!f$converged) {
        msg <- if (identical(f$failure, "sign change only across alpha = 0"))
          paste("the adjusted score has no regular root: its shape component changes sign only",
                "across alpha = 0, where the information is singular; the data do not support",
                "an estimate of the shape with this sample size")
        else paste("the adjusted score iteration did not converge:", f$failure)
        warning(msg, call. = FALSE)
      }
      theta <- f$theta
    }
    edf <- ncol(X) + 2
  }
  P <- length(theta)
  names(theta) <- c(colnames(Xf), "sigma", "alpha")
  V <- if (all(is.finite(theta))) {
    K <- sn_info(theta, Xf)
    H <- if (!is.null(S)) K + lam * S else K
    tryCatch(solve(H), error = function(e) matrix(NA_real_, P, P))
  } else matrix(NA_real_, P, P)
  dimnames(V) <- list(names(theta), names(theta))
  structure(list(coefficients = theta, vcov = V, method = method, call = cl,
                 X = Xf, y = y, p_lin = ncol(X), lambda = lam, edf = edf,
                 loglik = sn_loglik(theta, Xf, y), finite_mle = finite,
                 smooth = smooth_info, formula = formula),
            class = "snbr")
}

#' @export
print.snbr <- function(x, digits = 4, ...) {
  cat("Skew-normal regression, method:", x$method, "\n")
  if (!is.null(x$smooth)) cat("Smooth term in", x$smooth$var, " (edf of the model", format(x$edf, digits = 3), ")\n")
  th <- x$coefficients
  print(th[c(seq_len(x$p_lin), length(th) - 1L, length(th))], digits = digits)
  if (isFALSE(x$finite_mle)) cat("The maximum likelihood estimate of the shape is infinite.\n")
  invisible(x)
}

#' @export
summary.snbr <- function(object, ...) {
  th <- object$coefficients
  idx <- c(seq_len(object$p_lin), length(th) - 1L, length(th))
  se <- sqrt(diag(object$vcov))[idx]
  est <- th[idx]
  z <- est / se
  z[length(z) - 1L] <- NA   # no test for sigma = 0
  tab <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  structure(list(table = tab, object = object), class = "summary.snbr")
}

#' @export
print.summary.snbr <- function(x, ...) {
  print(x$object$call)
  cat("Method:", x$object$method, "\n\n")
  stats::printCoefmat(x$table, na.print = "")
  cat("\nLog-likelihood:", format(x$object$loglik, digits = 6),
      " (effective df", format(x$object$edf, digits = 3), ")\n")
  if (isFALSE(x$object$finite_mle)) cat("The maximum likelihood estimate of the shape is infinite.\n")
  invisible(x)
}

#' @export
coef.snbr <- function(object, ...) object$coefficients

#' @export
vcov.snbr <- function(object, ...) object$vcov

#' @export
logLik.snbr <- function(object, ...) {
  structure(object$loglik, df = object$edf, nobs = length(object$y), class = "logLik")
}

#' Skew-normal distribution function
#'
#' @param q quantiles.
#' @param xi,omega,alpha location, scale and shape.
#' @return Probabilities.
#' @export
psnorm <- function(q, xi = 0, omega = 1, alpha = 0) {
  z <- (q - xi) / omega
  n <- max(length(z), length(alpha))
  z <- rep_len(z, n)
  alpha <- rep_len(alpha, n)
  vapply(seq_len(n), function(i) {
    if (!is.finite(z[i])) return(as.numeric(z[i] > 0))
    a <- alpha[i]
    f <- function(t) 2 * stats::dnorm(t) * stats::pnorm(a * t)
    min(1, max(0, stats::integrate(f, -Inf, z[i], rel.tol = 1e-10)$value))
  }, numeric(1))
}

#' @export
residuals.snbr <- function(object, type = c("quantile", "response"), ...) {
  type <- match.arg(type)
  th <- object$coefficients
  P <- length(th)
  mu <- drop(object$X %*% th[seq_len(P - 2L)])
  if (type == "response") return(object$y - mu)
  stats::qnorm(psnorm(object$y, mu, th[P - 1L], th[P]))
}

#' @export
fitted.snbr <- function(object, ...) {
  th <- object$coefficients
  P <- length(th)
  drop(object$X %*% th[seq_len(P - 2L)])
}

#' @export
confint.snbr <- function(object, parm, level = 0.95, ...) {
  th <- object$coefficients
  se <- sqrt(diag(object$vcov))
  q <- stats::qnorm(1 - (1 - level) / 2)
  ci <- cbind(th - q * se, th + q * se)
  colnames(ci) <- paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2)), "%")
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
