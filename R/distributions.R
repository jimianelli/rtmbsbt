
#' Student's t-density
#' @param x A numeric value.
#' @param df Degrees of freedom.
#' @param mu Mean of the distribution.
#' @param sd Standard deviation of the distribution.
#' @param give_log Logical; if TRUE, returns the log density.
#' @return The density or log density of the Student's t-distribution.
dstudent <- function(x, df, mu, sd, give_log = FALSE) {
  sx <- (x - mu) / sd
  logres <- dt(sx, df, log = TRUE) - log(sd)
  if (give_log) logres else exp(logres)
}

#' Dirichlet density
#' @param x A numeric vector.
#' @param alpha A numeric vector of parameters.
#' @param give_log Logical; if TRUE, returns the log density.
#' @return The density or log density of the Dirichlet distribution.
ddirichlet <- function(x, alpha, give_log = FALSE) {
  logres <- lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
  if (give_log) logres else exp(logres)
}

#' Dirichlet-multinomial density
#' @param x A numeric vector.
#' @param alpha A numeric vector of parameters.
#' @param give_log Logical; if TRUE, returns the log density.
#' @return The density or log density of the Dirichlet-multinomial distribution.
ddm <- function(x, alpha, give_log = FALSE) {
  sum_alpha <- sum(alpha)
  logres <- lgamma(sum_alpha) - lgamma(sum(x) + sum_alpha) +
    sum(lgamma(x + alpha)) - sum(lgamma(alpha))
  if (give_log) logres else exp(logres)
}
