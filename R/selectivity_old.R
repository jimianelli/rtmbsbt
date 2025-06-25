
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
get_selectivity2 <- function(sel_switch_f, n_year, n_age, max_age, age_a, first_yr, first_yr_catch, 
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
