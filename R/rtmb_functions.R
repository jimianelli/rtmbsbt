# R/rtmb_functions.R

#' Transformation to ensure correlation parameters are between -1 and 1
#' @param x A numeric value.
#' @return A transformed value between -1 and 1.
rho_trans <- function(x) {
  2 / (1 + exp(-2 * x)) - 1
}

#' Inverse logit transformation with bounds
#' @param x A numeric value.
#' @param lb Lower bound (default is 0).
#' @param ub Upper bound (default is 1).
#' @return A transformed value between lb and ub.
invlogit_bnd <- function(x, lb = 0, ub = 1) {
  expx <- exp(x)
  (expx + lb) / (ub + expx)
}

#' Divide a vector by its mean
#' @param x A numeric vector.
#' @return A vector where each element is divided by the mean of x.
div_by_mean <- function(x) {
  x / mean(x)
}

#' Double logistic function
#' @param x A numeric vector.
#' @param p1 Parameter 1.
#' @param p2 Parameter 2.
#' @param p3 Parameter 3.
#' @return A transformed vector.
double_logistic <- function(x, p1, p2, p3) {
  gamma1 <- p1 + p2
  gamma2 <- p1 + gamma1 + p3
  (1 + exp(-log(19) * (x - gamma1) / p1))^(-1) *
    (1 - (1 + exp(-log(19) * (x - gamma2) / p3))^(-1)) *
    (0.95)^(-2)
}

#' Square of a value
#' @param x A numeric value.
#' @return The square of x.
square <- function(x) {
  x * x
}

#' Norm squared of a vector
#' @param x A numeric vector.
#' @return The sum of squares of the elements of x.
norm2 <- function(x) {
  sum(x * x)
}

#' First difference of a vector
#' @param x A numeric vector.
#' @return A vector representing the first difference of x.
first_difference <- function(x) {
  diff(x)
}

#' Third difference of a vector
#' @param x A numeric vector.
#' @return A vector representing the third difference of x.
third_difference <- function(x) {
  diff(diff(diff(x)))
}

#' Exponential minus one
#' @param x A numeric value.
#' @return exp(x) - 1, with a Taylor approximation for small x.
expm1 <- function(x) {
  if (abs(x) < 1e-5) {
    x + 0.5 * x * x
  } else {
    exp(x) - 1
  }
}

#' Posfun function
#' @param x A numeric value.
#' @param eps A small positive value.
#' @param pen A penalty term.
#' @return A value adjusted by posfun.
posfun <- function(x, eps, pen) {
  pen <- pen + ifelse(x < eps, 0.01 * (x - eps)^2, 0)
  ifelse(x >= eps, x, eps / (2 - x / eps))
}

#' Log-normal density
#' @param x A numeric value.
#' @param meanlog Mean of the log-normal distribution.
#' @param sdlog Standard deviation of the log-normal distribution.
#' @param give_log Logical; if TRUE, returns the log density.
#' @return The density or log density of the log-normal distribution.
dlnorm <- function(x, meanlog, sdlog, give_log = FALSE) {
  resid <- (log(x) - meanlog) / sdlog
  lp <- - (log(sqrt(2 * pi)) + 0.5 * resid^2 + log(x * sdlog))
  if (give_log) lp else exp(lp)
}

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

# R/rtmb_functions.R

#' Get selectivity skeleton
#' @param age_a A numeric vector of ages.
#' @param par_log_sel_skel A matrix of parameters for the selectivity skeleton.
#' @return An array of selectivity skeletons.
get_sel_skeleton <- function(age_a, par_log_sel_skel) {
  n_fishery <- 6
  sel_skeleton_fa <- array(list(), dim = c(n_fishery))
  for (f in 1:n_fishery) {
    sel_skeleton_fa[[f]] <- double_logistic(age_a, exp(par_log_sel_skel[f, 1]), 
                                           exp(par_log_sel_skel[f, 2]), exp(par_log_sel_skel[f, 3]))
  }
  return(sel_skeleton_fa)
}

#' Get selectivity
#' @param sel_switch_f A numeric vector indicating the selectivity switch for each fishery.
#' @param n_year The number of years.
#' @param n_age The number of ages.
#' @param max_age The maximum age.
#' @param age_a A numeric vector of ages.
#' @param first_yr The first year.
#' @param first_yr_catch The first year of catch.
#' @param first_yr_catch_f A numeric vector of the first year of catch for each fishery.
#' @param sel_min_age_f A numeric vector of the minimum age for selectivity for each fishery.
#' @param sel_max_age_f A numeric vector of the maximum age for selectivity for each fishery.
#' @param sel_end_f A numeric vector indicating if selectivity ends at the maximum age for each fishery.
#' @param sel_change_year_fy A matrix indicating the years when selectivity changes for each fishery.
#' @param sel_skeleton_fa An array of selectivity skeletons.
#' @param par_log_sel1_ay An array of parameters for selectivity 1.
#' @param par_log_sel2_ay An array of parameters for selectivity 2.
#' @param par_log_sel3_ay An array of parameters for selectivity 3.
#' @param par_log_sel4_ay An array of parameters for selectivity 4.
#' @param par_log_sel5_ay An array of parameters for selectivity 5.
#' @param par_log_sel6_ay An array of parameters for selectivity 6.
#' @return An array of selectivities.
get_selectivity <- function(sel_switch_f, n_year, n_age, max_age, age_a, first_yr, first_yr_catch, 
                            first_yr_catch_f, sel_min_age_f, sel_max_age_f, sel_end_f, 
                            sel_change_year_fy, sel_skeleton_fa, par_log_sel1_ay, par_log_sel2_ay, 
                            par_log_sel3_ay, par_log_sel4_ay, par_log_sel5_ay, par_log_sel6_ay) {
  
  n_fishery <- 6
  sel_zero <- rep(0, n_age)
  sel_fya <- array(list(), dim = c(n_fishery, n_year))
  
  for (f in 1:n_fishery) {
    for (y in 1:n_year) {
      sel_fya[[f, y]] <- sel_zero
    }
  }
  
  for (f in 1:n_fishery) {
    amin <- sel_min_age_f[f]
    amax <- sel_max_age_f[f]
    sel_tmp <- rep(0, amax - amin + 1)
    
    if (sel_switch_f[f] == 2) {
      ymin <- first_yr_catch - first_yr
      for (a in amin:amax) {
        if (f == 1) sel_tmp[a - amin + 1] <- exp(par_log_sel1_ay[a - amin + 1, 1])
        if (f == 2) sel_tmp[a - amin + 1] <- exp(par_log_sel2_ay[a - amin + 1, 1])
        if (f == 3) sel_tmp[a - amin + 1] <- exp(par_log_sel3_ay[a - amin + 1, 1])
        if (f == 4) sel_tmp[a - amin + 1] <- exp(par_log_sel4_ay[a - amin + 1, 1])
        if (f == 5) sel_tmp[a - amin + 1] <- exp(par_log_sel5_ay[a - amin + 1, 1])
        if (f == 6) sel_tmp[a - amin + 1] <- exp(par_log_sel6_ay[a - amin + 1, 1])
      }
      mean_tmp <- mean(sel_tmp)
      for (a in amin:amax) {
        sel_fya[[f, ymin]][a] <- sel_tmp[a - amin + 1] / mean_tmp
      }
      if (sel_end_f[f] && amax < max_age) {
        for (a in (amax + 1):n_age) {
          sel_fya[[f, ymin]][a] <- sel_fya[[f, ymin]][amax]
        }
      }
      ypar <- 2
      for (y in (ymin + 1):n_year) {
        if (sel_change_year_fy[f, y]) {
          for (a in amin:amax) {
            if (f == 1) sel_tmp[a - amin + 1] <- sel_fya[[f, y - 1]][a] * exp(par_log_sel1_ay[a - amin + 1, ypar])
            if (f == 2) sel_tmp[a - amin + 1] <- sel_fya[[f, y - 1]][a] * exp(par_log_sel2_ay[a - amin + 1, ypar])
            if (f == 3) sel_tmp[a - amin + 1] <- sel_fya[[f, y - 1]][a] * exp(par_log_sel3_ay[a - amin + 1, ypar])
            if (f == 4) sel_tmp[a - amin + 1] <- sel_fya[[f, y - 1]][a] * exp(par_log_sel4_ay[a - amin + 1, ypar])
            if (f == 5) sel_tmp[a - amin + 1] <- sel_fya[[f, y - 1]][a] * exp(par_log_sel5_ay[a - amin + 1, ypar])
            if (f == 6) sel_tmp[a - amin + 1] <- sel_fya[[f, y - 1]][a] * exp(par_log_sel6_ay[a - amin + 1, ypar])
          }
          ypar <- ypar + 1
          mean_tmp <- mean(sel_tmp)
          for (a in amin:amax) {
            sel_fya[[f, y]][a] <- sel_tmp[a - amin + 1] / mean_tmp
          }
          if (sel_end_f[f] && amax < max_age) {
            for (a in (amax + 1):n_age) {
              sel_fya[[f, y]][a] <- sel_fya[[f, y]][amax]
            }
          }
        } else {
          sel_fya[[f, y]] <- sel_fya[[f, y - 1]]
        }
      }
    } else {
      ymin <- first_yr_catch_f[f] - first_yr
      ypar <- 1
      for (y in ymin:n_year) {
        if (sel_change_year_fy[f, y]) {
          for (a in amin:amax) {
            if (f == 1) sel_tmp[a - amin + 1] <- sel_skeleton_fa[[f]][a] + par_log_sel1_ay[a - amin + 1, ypar]
            if (f == 2) sel_tmp[a - amin + 1] <- sel_skeleton_fa[[f]][a] + par_log_sel2_ay[a - amin + 1, ypar]
            if (f == 3) sel_tmp[a - amin + 1] <- sel_skeleton_fa[[f]][a] + par_log_sel3_ay[a - amin + 1, ypar]
            if (f == 4) sel_tmp[a - amin + 1] <- sel_skeleton_fa[[f]][a] + par_log_sel4_ay[a - amin + 1, ypar]
            if (f == 5) sel_tmp[a - amin + 1] <- sel_skeleton_fa[[f]][a] + par_log_sel5_ay[a - amin + 1, ypar]
            if (f == 6) sel_tmp[a - amin + 1] <- sel_skeleton_fa[[f]][a] + par_log_sel6_ay[a - amin + 1, ypar]
          }
          ypar <- ypar + 1
          sel_tmp <- exp(sel_tmp)
          mean_tmp <- mean(sel_tmp)
          for (a in amin:amax) {
            sel_fya[[f, y]][a] <- sel_tmp[a - amin + 1] / mean_tmp
          }
          if (sel_end_f[f] && amax < max_age) {
            for (a in (amax + 1):n_age) {
              sel_fya[[f, y]][a] <- sel_fya[[f, y]][amax]
            }
          }
        } else {
          sel_fya[[f, y]] <- sel_fya[[f, y - 1]]
        }
      }
    }
  }
  
  return(sel_fya)
}

#' Construct precision matrix Q
#' @param n_years The number of years.
#' @param n_ages The number of ages.
#' @param ay_Index A matrix of indices.
#' @param rho_y Partial correlation by years.
#' @param rho_a Partial correlation by ages.
#' @param rho_c Partial correlation by cohort.
#' @param log_sigma2 Variance parameter governing GMRF.
#' @param Var_Param Parameterization of variance (0=conditional, 1=marginal).
#' @return A sparse precision matrix.
construct_Q <- function(n_years, n_ages, ay_Index, rho_y, rho_a, rho_c, log_sigma2, Var_Param) {
  total_n <- n_years * n_ages
  B <- Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(total_n, total_n))
  I <- Matrix::sparseMatrix(i = 1:total_n, j = 1:total_n, x = 1, dims = c(total_n, total_n))
  Omega <- Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(total_n, total_n))
  Q_sparse <- Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(total_n, total_n))
  
  for (n in 1:total_n) {
    age <- ay_Index[n, 1]
    year <- ay_Index[n, 2]
    if (age > 1) {
      for (n1 in 1:total_n) {
        if (ay_Index[n1, 1] == age - 1 && ay_Index[n1, 2] == year) {
          B[n, n1] <- rho_y
        }
      }
    }
    if (year > 1) {
      for (n1 in 1:total_n) {
        if (ay_Index[n1, 1] == age && ay_Index[n1, 2] == year - 1) {
          B[n, n1] <- rho_a
        }
      }
    }
    if (year > 1 && age > 1) {
      for (n1 in 1:total_n) {
        if (ay_Index[n1, 1] == age - 1 && ay_Index[n1, 2] == year - 1) {
          B[n, n1] <- rho_c
        }
      }
    }
  }
  
  if (Var_Param == 0) {
    for (i in 1:total_n) {
      Omega[i, i] <- 1 / exp(log_sigma2)
    }
  }
  
  if (Var_Param == 1) {
    L <- solve(I - B)
    d <- rep(0, total_n)
    for (n in 1:total_n) {
      if (n == 1) {
        d[n] <- exp(log_sigma2)
      } else {
        cumvar <- 0
        for (n1 in 1:(n - 1)) {
          cumvar <- cumvar + L[n, n1] * d[n1] * L[n, n1]
        }
        d[n] <- (exp(log_sigma2) - cumvar) / (L[n, n]^2)
      }
    }
    for (i in 1:total_n) {
      Omega[i, i] <- 1 / d[i]
    }
  }
  
  B_transpose <- Matrix::t(B)
  Q_sparse <- (I - B_transpose) %*% Omega %*% (I - B)
  return(Q_sparse)
}
