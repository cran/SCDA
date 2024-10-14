#' @keywords internal
#' @noRd

########## semlogLik_i: computes the log-likelihood function for the SEM model
semlogLik_i <- function(Y, X, W, beta, lamb, sigma2) {
  N <- length(Y)  # Number of observations
  I <- diag(N)    # Identity matrix
  W <- as.matrix(W)

  # Compute the log determinant term
  log_det_term <- log(det(as.matrix((I - lamb * W))))

  # Compute the residuals
  residuals <- as.matrix(I - lamb * W) %*% (Y - X %*% beta)

  # Compute the log-likelihood
  # Errore
  # log_likelihood <- - (1 / 2) * log(2 * pi * sigma2) + (log_det_term/N) - (1 / (2 * sigma2)) * residuals^2
  # Corretta (3.14 di LeSage-Pace 2009)
  # log_likelihood <- - (N / 2) * log(pi * sigma2) + log_det_term - (residuals^2 / (2 * sigma2))
  # Average log-likelihood: provides the average contribution of the single observation and numerical stability
  log_likelihood <- - (1 / 2) * log(2 * pi * sigma2) + log_det_term/N - (residuals^2 / (2 * sigma2))

  return(log_likelihood)
}
