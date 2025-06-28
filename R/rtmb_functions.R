#' Helper to make closure
#' 
#' @param f a \code{vector} of midpoints.
#' @param d natural mortality at the reference size.
#' @return a \code{vector}.
#' @export
#' 
cmb <- function(f, d) function(p) f(p, d)


#' Tag Recapture Likelihood
#'
#' Computes the negative log-likelihood for tag-recapture data using a seasonal, age-based model.
#'
#' @param tag_switch Integer flag to turn the likelihood on (>0) or off (0).
#' @param minK Integer, starting cohort index.
#' @param n_K Number of cohorts.
#' @param n_T Number of tagging time periods.
#' @param n_I Number of release ages.
#' @param n_J Number of recapture ages.
#' @param first_yr First model year.
#' @param M_a Vector of natural mortality at age.
#' @param hrate_ysa 3D array [year, season, age] of harvest rates.
#' @param par_hstar_i Vector of free parameters for tag shedding (h*) for each cohort and release age.
#' @param tag_release_cta 3D array [cohort, time, age] of numbers of tags released.
#' @param tag_recap_ctaa 4D array [cohort, time, rel_age, recap_age] of numbers of tags recaptured.
#' @param minI,maxI,maxJ Vectors giving the min and max release and recapture ages for each cohort.
#' @param shed1,shed2 Vectors of shedding rates per time period.
#' @param tag_rep_rates_ya Matrix of tag reporting rates by year and age.
#' @param tag_H_factor Numeric scaling factor for incomplete mixing.
#' @param tag_var_factor Overdispersion factor.
#' @return Negative log-likelihood (scalar) for tag data.
#' @export
#'
get_tag_like <- function(tag_switch, minK, n_K, n_T, n_I, n_J, first_yr, M_a, hrate_ysa, 
                         tag_release_cta, tag_recap_ctaa, 
                         minI, maxI, maxJ, shed1, shed2, tag_rep_rates_ya,
                         tag_H_factor, tag_var_factor) {
  
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  prR <- tag_pred <- tag_resid <- array(0, dim = c(n_K, n_T, n_I, n_J))
  S1 <- S2 <- f1 <- f2 <- array(0, dim = c(n_K, n_T, n_J))
  S1star <- S2star <- tag_release_adj <- array(0, dim = c(n_K, n_T, n_I))
  
  ## calculate number of independent Dirichlet-multinomial likelihood components, 
  ## i.e., one for every cohort, tagger group and release age, but leaving out those with zero releases 
  ## and for which the max recapture age = release age (otherwise it would just be n_K * n_T * n_I)
  ntaglike <- 0
  for (k in seq_len(n_K)) {
    for (t in seq_len(n_T)) {
      for (i in minI[k]:maxI[k]) {
        if ((tag_release_cta[k, t, i] > 0) & (maxJ[k] > i)) ntaglike <- ntaglike + 1
      }
    }
  }
  
  # Survival and exploitation rates for tagged fish by age and tagger
  for (k in seq_len(n_K)) {
    for (j in minI[k]:maxJ[k]) {
      iy <- minK + k - 1 + j
      hplus <- tag_H_factor * hrate_ysa[iy, 1, j + 1] # factor to account for lack of complete mixing in season 1 
      for (t in seq_len(n_T)) {
        S1[k, t, j] <- (1 - hplus) * (1 - hrate_ysa[iy, 2, j + 1]) * exp(-M_a[j + 1] - shed2[t])
        S2[k, t, j] <- (1 - hplus) * (1 - hrate_ysa[iy, 2, j + 1]) * exp(-M_a[j + 1] - 2 * shed2[t])
        f1[k, t, j] <- hplus + (1 - hplus) * exp(-0.5 * (M_a[j + 1] + shed2[t])) * hrate_ysa[iy, 2, j + 1]
        f2[k, t, j] <- hplus + (1 - hplus) * exp(-0.5 * (M_a[j + 1] + 2 * shed2[t])) * hrate_ysa[iy, 2, j + 1]
      }
    }
  }
  
  for (k in seq_len(n_K)) {
    for (i in minI[k]:maxI[k]) {
      for (t in seq_len(n_T)) {
        iy <- minK + k - 1 + i
        S1star[k, t, i] <-  exp(-M_a[i + 1] - shed2[t])
        S2star[k, t, i] <-  exp(-M_a[i + 1] - 2 * shed2[t])
        # calculate adjusted tag release numbers by subtracting number of recaptures in year of tagging (taking into account non-reporting)
        tag_release_adj[k, t, i] <- tag_release_cta[k, t, i] - tag_recap_ctaa[k, t, i, i] / tag_rep_rates_ya[k + i - 2, i]
      }
    }
  }
  
  # Generate probabilities of recapture
  # tag_rep_rates_ya: years 1991-1997, ages 1-8
  for (k in seq_len(n_K)) {
    for (t in seq_len(n_T)) {
      for (i in minI[k]:maxI[k]) {
        for (j in i:maxJ[k]) {
          # no longer include recaptures in year of tagging in likelihood  so deleted case where j = i
          if (j == (i + 1)) {
            prR[k, t, i, j] <- (2 * shed1[t] * S1star[k, t, i] * f1[k, t, j] - shed1[t]^2 * S2star[k, t, i] * f2[k, t, j]) * tag_rep_rates_ya[k + j - 2, j]
          }
          if (j > (i + 1)) {
            prodS1 <- prodS2 <- 1
            for (s in (i + 1):(j - 1)) {
              prodS1 <- prodS1 * S1[k, t, s]
              prodS2 <- prodS2 * S2[k, t, s]
            }
            prR[k, t, i, j] <- (2 * shed1[t] * S1star[k, t, i] * prodS1 * f1[k, t, j] - 
                                  shed1[t]^2 * S2star[k, t, i] * prodS2 * f2[k, t, j]) * tag_rep_rates_ya[k + j - 2, j]
          }
        }
      }
    }
  }
  
  # Predicted numbers of recaptures and residuals
  for (k in seq_len(n_K)) {
    for (t in seq_len(n_T)) {
      for (i in minI[k]:maxI[k]) {
        for (j in (i + 1):maxJ[k]) {
          tag_pred[k, t, i, j] <- tag_release_adj[k, t, i] * prR[k, t, i, j] 
          # THE TAG RESID IS WRONG AND NEEDS TO BE OSA RESIDS
          tag_resid[k, t, i, j] <- (tag_recap_ctaa[k, t, i, j] - tag_pred[k, t, i, j]) / 
            sqrt(1e-5 + tag_var_factor * tag_pred[k, t, i, j] * (1 - prR[k, t, i, j]))
        }
      }
    }
  }
  
  # 8) log-likelihood under beta-binomial
  lp <- numeric(ntaglike)
  if (tag_switch > 0) {
    index <- 1
    for (k in seq_len(n_K)) {
      for (t in seq_len(n_T)) {
        for (i in minI[k]:maxI[k]) {
          if (tag_release_cta[k,t,i] > 0 & maxJ[k] > i) {
            loglkhd_R <- 0 # loglkhd for cohort k, tagger t and release age i (summed over corresponding recaptures )
            tag_od <- (tag_release_cta[k,t,i] - tag_var_factor) / (tag_var_factor - 1)
            if (tag_od < 0) tag_od <- 1e-3
            loglkhd_R <- loglkhd_R + lgamma(tag_od) - lgamma(tag_release_adj[k,t,i] + tag_od)
            totR <- totprR <- 0
            for (j in (i + 1):maxJ[k]) {
              loglkhd_R <- loglkhd_R + lgamma(tag_recap_ctaa[k,t,i,j] + tag_od * prR[k,t,i,j]) - lgamma(tag_od * prR[k,t,i,j]) 
              totR <- totR + tag_recap_ctaa[k,t,i,j]
              totprR <- totprR + prR[k,t,i,j]
            }
            notR <- tag_release_adj[k,t,i] - totR
            pr_notR <- 1 - totprR
            loglkhd_R <- loglkhd_R + lgamma(notR + tag_od * pr_notR) - lgamma(tag_od * pr_notR)
            lp[index] <- -loglkhd_R
            index <- index + 1
          }
        }
      }
    }
  }
  return(list(pred = tag_pred, resid = tag_resid, lp = lp))
}

#' Aerial Survey Likelihood
#'
#' Computes the likelihood contribution from an aerial index of abundance survey.
#'
#' @param aerial_switch Integer specifying the type of selectivity (0–4).
#' @param aerial_years Vector of years corresponding to observations.
#' @param aerial_obs Vector of observed aerial survey indices.
#' @param aerial_cv Coefficients of variation for observations.
#' @param aerial_cov Covariance matrix of aerial observations.
#' @param par_aerial_tau Observation error term.
#' @param par_log_aerial_sel Vector of log selectivity parameters (for ages 2 and 4).
#' @param number_ysa 3D array [year, season, age] of numbers-at-age.
#' @param weight_fya 3D array [fleet, year, age] of weights-at-age.
#' @return A list with predicted values, residuals, and log-likelihood components:
#' \describe{
#'   \item{pred}{Predicted survey index.}
#'   \item{resid}{Log residuals.}
#'   \item{lp_aerial_tau}{Log determinant of the covariance matrix.}
#'   \item{lp}{Negative log-likelihood contribution.}
#' }
#' @importFrom RTMB ADoverload dmvnorm
#' @export
#'
get_aerial_survey_like <- function(aerial_switch, aerial_years, aerial_obs, aerial_cv, aerial_cov, 
                                   par_aerial_tau, par_log_aerial_sel, number_ysa, weight_fya) {
  
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_aerial <- length(aerial_obs)
  aerial_tau2 <- par_aerial_tau^2
  I <- diag(nrow = n_aerial, ncol = n_aerial)
  cov_matrix <- aerial_cov + I * aerial_tau2
  cov_inv <- solve(cov_matrix)
  
  # Aerial selectivity (ages 2, 3, and 4)
  if (aerial_switch == 0) {
    aerial_sel <- c(0.5, 1, 1)
  } else if (aerial_switch == 1) {
    aerial_sel <- c(exp(par_log_aerial_sel[1]), 1, exp(par_log_aerial_sel[2]))
  } else if (aerial_switch == 2) {
    aerial_sel <- rep(1, 3)
  } else if (aerial_switch == 3) {
    aerial_sel <- c(0.33, 1, 0.33)
  } else if (aerial_switch == 4) {
    aerial_sel <- c(0.5, 1, 1)
  } else {
    stop("Unknown aerial_switch")
  }
  aerial_pred <- numeric(n_aerial)
  axx <- 1:3
  for (i in seq_len(n_aerial)) {
    y <- aerial_years[i]
    aerial_pred[i] <- sum(number_ysa[y, 1, axx + 2] * weight_fya[6, y, axx + 2] * aerial_sel[axx])
    #for (a in 1:3) { # Ages 2,3,4 = indices 3,4,5 (1-based)
    #  aerial_pred[i] <- aerial_pred[i] + number_ysa[y, 1, a + 2] * weight_fya[6, y, a + 2] * aerial_sel[a]
    #}
  }
  aerial_resid <- log(aerial_obs) - log(aerial_pred) # Residuals
  aerial_log_q <- sum(cov_inv %*% aerial_resid) / sum(cov_inv) # Estimate log_q
  aerial_resid <- aerial_resid - aerial_log_q
  aerial_pred <- aerial_pred * exp(aerial_log_q) # Scale predictions
  aerial_log_pred <- log(aerial_pred)
  aerial_log_obs <- log(aerial_obs)
  lp_aerial_tau <- 0.5 * log(det(cov_matrix))
  if (aerial_switch > 0) {
    # lp <- 0.5 * as.vector(aerial_resid %*% cov_inv %*% aerial_resid) # Mahalanobis term: 0.5 * x' Sigma^{-1} x
    lp <- -dmvnorm(aerial_log_obs, aerial_log_pred, cov_matrix, log = TRUE)
  } else {
    lp <- 0
  }
  return(list(pred = aerial_pred, resid = aerial_resid, lp_aerial_tau = lp_aerial_tau, lp = lp))
}


#' Parent-Offspring Pairing Likelihood
#'
#' Computes the negative log-likelihood for observed parent-offspring pairings from genetic samples.
#'
#' @param pop_switch Integer flag for activation (currently unused).
#' @param pop_obs Matrix [n,6]: capture year, offspring year, adult age/length, type flag (0=age, 1=length), number of pairings, and number of comparisons.
#' @param phi_ya Matrix [year, age] of reproductive output-at-age.
#' @param paly Array [length, age, year] of predicted adult distributions at length.
#' @param spawning_biomass_y Vector of spawning biomass by year.
#' @return Vector of negative log-likelihood contributions for each observation.
#' @importFrom RTMB ADoverload
#' @export
#'
get_POP_like <- function(pop_switch, pop_obs, phi_ya, paly, spawning_biomass_y) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_pops <- nrow(pop_obs)
  n_age <- ncol(phi_ya)
  lp <- numeric(n_pops)
  for (i in seq_len(n_pops)) {
    cc <- pop_obs[i, 1]
    yy <- pop_obs[i, 2]
    nP <- pop_obs[i, 5]
    nC <- pop_obs[i, 6]
    if (pop_obs[i, 4] == 0) { # Adult has a direct age
      aa <- pop_obs[i, 3] + 1 # Estimated age of adult at time of capture
      ba <- ifelse(aa - (yy - cc) < 1, 1, aa - (yy - cc))
      ba <- min(ba, n_age)
      pp <- (2 * phi_ya[cc, ba]) / spawning_biomass_y[cc] # Parental probability
    }
    if (pop_obs[i, 4] == 1) { # Adult has observed length only
      ll <- pop_obs[i, 3] # Observed length bin of adult at time of capture
      amin <- yy - cc + 1 # Anything younger than this can't be a parent
      arng <- amin:n_age
      ba <- arng - (yy - cc)
      pp <- (2 / spawning_biomass_y[cc]) * sum(paly[ll, arng, yy] * phi_ya[cc, ba])
    }
    lp[i] <- -(nP * log(pp) + (nC - nP) * log(1 - pp))
  }
  return(lp)
}

#' Half-Sibling Pair Likelihood
#'
#' Computes the negative log-likelihood for half-sibling pair observations based on reproduction and survival.
#'
#' @param hsp_switch Integer flag for activation (currently unused).
#' @param hsp_obs Matrix [n,5] with columns: cohort1, cohort2, cohort difference, comparisons, matches.
#' @param hsp_q Effective sampling fraction of potential parents.
#' @param hsp_false_negative Probability of detecting a true match (1 – false negative rate).
#' @param number_ysa 3D array [year, season, age] of numbers-at-age.
#' @param phi_ya Matrix [year, age] of reproductive output-at-age.
#' @param M_a Vector of natural mortality at age.
#' @param spawning_biomass_y Vector of spawning biomass by year.
#' @param hrate_ysa 3D array [year, season, age] of harvest rates.
#' @return Vector of negative log-likelihood contributions for each observation.
#' @importFrom RTMB ADoverload dbinom
#' @export
#'
get_HSP_like <- function(hsp_switch, hsp_obs, hsp_q, hsp_false_negative, 
                         number_ysa, phi_ya, M_a, spawning_biomass_y, hrate_ysa) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_hsp <- nrow(hsp_obs)
  n_year <- nrow(phi_ya)
  n_age <- ncol(phi_ya)
  lp <- numeric(n_hsp)
  nK <- hsp_obs[, 5]
  # nK <- OBS(nK) # tag it for simulation
  # Calculate the relative age distribution in the adults for each year and age
  gamx_ya <- matrix(0, nrow = n_year, ncol = n_age)
  for (y in seq_len(n_year)) {
    gamx_ya[y,] <- number_ysa[y, 2,] * phi_ya[y,] / spawning_biomass_y[y] # Use season 2 for numbers-at-age of adults
  }
  for (i in seq_len(n_hsp)) {
    cmin <- hsp_obs[i, 1] + 1
    cmax <- hsp_obs[i, 2] + 1
    cdif <- hsp_obs[i, 3]
    nC <- hsp_obs[i, 4]
    xtmp <- 0
    for (a in seq_len(n_age)) {
      # Step 2: calculate cumulative survival between c1 and c2 given reference adult age
      cumS <- 1
      if ((a + cdif) < n_age) {
        for (ia in seq_len(cdif)) {
          age_idx <- a + ia - 1
          year_idx <- cmin + ia - 1
          cumS <- cumS * exp(-M_a[age_idx]) * (1 - hrate_ysa[year_idx, 1, age_idx]) * (1 - hrate_ysa[year_idx, 2, age_idx])
        }
      } else {
        for (ia in seq_len(cdif)) {
          age_idx <- a + ia - 1
          year_idx <- cmin + ia - 1
          idx <- ifelse(age_idx < n_age, age_idx, n_age)
          cumS <- cumS * exp(-M_a[idx]) * (1 - hrate_ysa[year_idx, 1, idx]) * (1 - hrate_ysa[year_idx, 2, idx])
        }
      }
      # Step 3: effective increase in RO-at-age
      if ((a + cdif) <= n_age) {
        phi_val <- phi_ya[cmax, a + cdif]
      } else {
        phi_val <- phi_ya[cmax, n_age]
      }
      # Step 4: integration
      xtmp <- xtmp + gamx_ya[cmin, a] * cumS * phi_val
    }
    pp <- 4 * hsp_q * xtmp / spawning_biomass_y[cmax]
    phsp <- pp * hsp_false_negative
    if (hsp_switch > 0) lp[i] <- -dbinom(nK[i], nC, pp, log = TRUE)
    # if (hsp_switch > 0 && phsp > 0) { if (nK > 0) {lp[i] <- -(nK * log(phsp) + (nC - nK) * log(1 - phsp)) } else lp[i] <- -nC * log(1 - phsp)
  }
  return(lp)
}


#' Selectivity Prior Penalty Likelihood
#'
#' Computes penalty terms for initial selectivity levels, year-to-year changes, and smoothness constraints.
#'
#' @param first_yr Model start year.
#' @param first_yr_catch_f Vector of first years of catch data per fishery.
#' @param sel_min_age_f, sel_max_age_f Vectors of minimum and maximum selectivity ages per fishery.
#' @param sel_change_year_fy Logical matrix indicating whether selectivity changed in a year.
#' @param sel_change_sd_fy Matrix of standard deviations for selectivity changes.
#' @param sel_smooth_sd_f Vector of smoothing SDs per fishery.
#' @param par_sels_init_i Vector of initial selectivity log-values.
#' @param par_sels_change_i Vector of selectivity changes.
#' @param sel_fya 3D array [fishery, year, age] of selectivity values.
#' @return A numeric vector of 3 penalty components: change penalty, smoothness penalty, and log-mean prior.
#' @export
#'
get_sel_like <- function(first_yr, first_yr_catch_f, sel_min_age_f, sel_max_age_f, sel_change_year_fy, sel_change_sd_fy, 
                         sel_smooth_sd_f, par_sels_init_i, par_sels_change_i, sel_fya) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_fishery <- dim(sel_fya)[1]
  n_year <- dim(sel_fya)[2]
  n_age <- dim(sel_fya)[3]
  lp <- numeric(3)
  ipar <- 1
  jpar <- 1
  for (f in seq_len(n_fishery)) {
    amin <- sel_min_age_f[f] + 1
    amax <- sel_max_age_f[f] + 1
    n_sel <- amax - amin + 1
    # Initial selectivity prior (penalize log mean)
    sel_tmp <- exp(par_sels_init_i[ipar:(ipar + n_sel - 1)])
    lp[3] <- lp[3] + 50 * log(mean(sel_tmp))^2
    ipar <- ipar + n_sel
    for (y in seq(from = first_yr_catch_f[f] - first_yr + 1, to = n_year)) {
      # Selectivity change penalty
      if (sel_change_year_fy[f, y]) {
        sel_change <- par_sels_change_i[jpar:(jpar + n_sel - 1)] / sel_change_sd_fy[f, y]
        lp[1] <- lp[1] + 0.5 * sum(sel_change^2)
        jpar <- jpar + n_sel
      }
      # Smoothness penalty (third-difference on log-scale)
      smooth_sd <- ifelse(f == 3, 20, sel_smooth_sd_f[f])
      sel_log <- log(sel_fya[f, y, amin:amax])
      if (n_sel >= 4) {
        td <- diff(diff(diff(sel_log)))
        lp[2] <- lp[2] + 0.5 * sum((td / smooth_sd)^2)
      }
    }
  }
  return(lp)
}

#' Recruitment prior
#'
#' Applies a split variance penalty to recruitment deviations.
#'
#' @param rdev_y Vector of recruitment deviations.
#' @param sigma_r Recruitment standard deviation.
#' @param tau_ac2 Temporal autocorrelation squared.
#' @return Negative log-prior penalty (scalar).
#' @export
#'
get_recruitment_prior <- function(rdev_y, sigma_r, tau_ac2) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_year <- length(rdev_y)
  r1 <- rdev_y[1:(n_year - 3)]
  r2 <- rdev_y[(n_year - 2):n_year]
  lp <- n_year * log(sigma_r) + 
    0.5 * sum(r1^2) / sigma_r^2 + 
    0.5 * sum(r2^2) / (sigma_r^2 * (1 - tau_ac2^2))
  return(lp)
}


#' CPUE Index Likelihood
#'
#' Computes the likelihood for a standardized CPUE index using a log-linear model.
#'
#' @param cpue_switch Integer switch to activate the likelihood.
#' @param cpue_a1, cpue_a2 Minimum and maximum CPUE age indices.
#' @param cpue_years Vector of year indices for CPUE observations.
#' @param cpue_obs Observed CPUE values.
#' @param cpue_adjust Vector of adjustment scalars for each CPUE observation.
#' @param cpue_sigma Observation standard deviation.
#' @param cpue_omega Power parameter for scaling to total numbers.
#' @param log_cpue_q Logarithm of catchability coefficient.
#' @param number_ysa 3D array [year, season, age] of numbers-at-age.
#' @param sel_fya 3D array [fishery, year, age] of selectivity.
#' @return List with predicted CPUE, residuals, and likelihood vector.
#' @export
get_cpue_like <- function(cpue_switch, cpue_a1 = 5, cpue_a2 = 17, 
                          cpue_years, cpue_obs, cpue_adjust, cpue_sigma, cpue_omega, 
                          log_cpue_q, par_cpue_creep, number_ysa, sel_fya) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_cpue <- length(cpue_obs)
  n_age <- dim(sel_fya)[3]
  cpue_adjust <- cpue_log_pred <- numeric(n_cpue)
  cpue_adjust[1] <- 1
  for (i in 2:n_cpue) cpue_adjust[i] <- cpue_adjust[i - 1] + par_cpue_creep
  for (i in seq_len(n_cpue)) {
    y <- cpue_years[i]
    cpue_sel <- sel_fya[7, y, 5:n_age] # age 4+
    cpue_n <- number_ysa[y, 2, 5:n_age] # season 2
    cpue_selm <- sel_fya[7, y, (cpue_a1 + 1):(cpue_a2 + 1)]
    tmpN <- sum(cpue_sel * cpue_n) / mean(cpue_selm)
    cpue_log_pred[i] <- log(cpue_adjust[i]) + cpue_omega * log(tmpN)
  }
  cpue_log_pred <- cpue_log_pred - log(mean(exp(cpue_log_pred))) + log_cpue_q
  cpue_log_obs <- log(cpue_obs)
  cpue_pred <- exp(cpue_log_pred)
  cpue_resid <- (log(cpue_obs) - cpue_log_pred) / cpue_sigma
  if (cpue_switch > 0) {
    #lp <- log(cpue_sigma) + 0.5 * cpue_resid^2
    lp <- -dnorm(cpue_log_obs, cpue_log_pred, cpue_sigma, log = TRUE)
  } else {
    lp <- numeric(n_cpue)
  }
  return(list(pred = cpue_pred, resid = cpue_resid, lp = lp))
}

#' Age composition likelihood
#'
#' Calculates the multinomial likelihood for observed age compositions.
#'
#' @param removal_switch_f Vector of year indices.
#' @param af_year Vector of year indices.
#' @param af_fishery Vector of fishery indices.
#' @param af_min_age, af_max_age Vectors of minimum and maximum ages.
#' @param af_obs Matrix of observed proportions-at-age.
#' @param af_n Vector of effective sample sizes.
#' @param catch_pred_fya 3D array [fishery, year, age] of predicted catch.
#' @return List with predicted age compositions and negative log-likelihoods.
#' @export
#' 
get_age_like <- function(removal_switch_f, af_year, af_fishery, af_min_age, af_max_age, af_obs, af_n, catch_pred_fya) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_af <- nrow(af_obs)
  n_age <- dim(catch_pred_fya)[3]
  lp <- numeric(n_af)
  age_pred <- matrix(0, n_af, n_age)
  for (i in seq_len(n_af)) {
    f <- af_fishery[i]
    y <- af_year[i]
    amin <- af_min_age[i] + 1
    amax <- af_max_age[i] + 1
    n_a <- amax - amin + 1
    obs <- af_obs[i, amin:amax]
    obs <- (obs / sum(obs)) + 1e-6
    pred <- catch_pred_fya[f, y, amin:amax]
    pred[1] <- sum(catch_pred_fya[f, y, 1:amin])
    pred[n_a] <- sum(catch_pred_fya[f, y, amax:n_age])
    pred <- (pred / sum(pred)) + 1e-6
    age_pred[i, amin:amax] <- pred
    if (removal_switch_f[f] == 0) {
      lp[i] <- af_n[i] * sum(obs * log(obs)) - af_n[i] * sum(obs * log(pred))
    }
  }
  return(list(pred = age_pred, lp = lp))
}


#' Length Composition Likelihood
#'
#' Computes likelihood for observed length compositions using ALKs.
#'
#' @param removal_switch_f, lf_season, lf_fishery Vectors of indices for observations.
#' @param lf_year, lf_season, lf_fishery Vectors of indices for observations.
#' @param lf_minbin Minimum size bin to be aggregated.
#' @param lf_obs Matrix of observed length proportions.
#' @param lf_n Vector of effective sample sizes.
#' @param catch_pred_fya 3D array of predicted catch.
#' @param alk_ysal 4D array [year, season, age, length_bin] of ALKs.
#' @return List with predicted compositions and likelihood contributions.
#' @export
#' 
get_length_like <- function(removal_switch_f, lf_year, lf_season, lf_fishery, lf_minbin, lf_obs, 
                            lf_n, catch_pred_fya, alk_ysal) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_lf <- nrow(lf_obs)
  n_bins <- 25
  n_age <- dim(catch_pred_fya)[3]
  lp <- numeric(n_lf)
  lf_pred <- matrix(0, n_lf, n_bins)
  for (i in seq_len(n_lf)) {
    y <- lf_year[i]
    s <- lf_season[i]
    f <- lf_fishery[i]
    mbin <- lf_minbin[f]
    catch_a <- catch_pred_fya[f, y,]
    pred <- catch_a %*% alk_ysal[y, s,, 1:n_bins]
    pred <- pred[1,] / sum(pred)
    obs <- lf_obs[i,]
    if (mbin > 1) {
      pred[mbin] <- sum(pred[1:mbin])
      obs[mbin] <- sum(obs[1:mbin])
      obs <- obs[mbin:n_bins]
      pred <- pred[mbin:n_bins]
    }
    lf_pred[i, mbin:n_bins] <- pred
    if (removal_switch_f[f] == 0) {
      lp[i] <- -lf_n[i] * sum(obs * log(pred))
      obs <- obs + 1e-6
      lp[i] <- lp[i] + lf_n[i] * sum(obs * log(obs))
    }
  }
  return(list(pred = lf_pred, lp = lp))
}


get_cpue_length_like <- function(cpue_years, cpue_lfs, cpue_n, number_ysa, sel_fya, alk_ysal) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_lf <- nrow(cpue_lfs)
  n_age <- dim(sel_fya)[3]
  n_bins <- 25
  lp <- numeric(n_lf)
  lf_pred <- matrix(0, n_lf, n_bins)
  for (i in seq_len(n_lf)) {
    y <- cpue_years[i]
    catch_a <- number_ysa[y, 2,] * sel_fya[7, y,]
    pred <- catch_a %*% alk_ysal[y, 2,, 1:n_bins]
    pred <- pred[1,] / sum(pred)
    lf_pred[i,] <- pred
    obs <- cpue_lfs[i,]
    lp[i] <- -cpue_n[i] * sum(obs * log(pred))
    obs <- obs + 1e-6
    lp[i] <- lp[i] + cpue_n[i] * sum(obs * log(obs))
  }
  return(list(pred = lf_pred, lp = lp))
}

#' Troll CPUE Likelihood
#'
#' Computes the log-likelihood for troll CPUE data.
#'
#' @param troll_switch Integer flag to turn likelihood on or off.
#' @param troll_years Vector of year indices.
#' @param troll_obs Observed troll index values.
#' @param troll_sd Observation standard errors.
#' @param troll_tau Process error.
#' @param number_ysa 3D array of numbers-at-age.
#' @return List of predicted values, residuals, and log-likelihoods.
#' @export
#'
get_troll_like <- function(troll_switch, troll_years, troll_obs, troll_sd, troll_tau, number_ysa) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_troll <- length(troll_years)
  troll_sd2 <- troll_sd^2 + troll_tau^2
  pred <- lp <- numeric(n_troll)
  for (i in seq_len(n_troll)) {
    iy <- troll_years[i]
    pred[i] <- number_ysa[iy, 1, 2]
  }
  troll_q <- sum(pred * (troll_obs / troll_sd2)) / sum(pred * (pred / troll_sd2))
  troll_pred <- pred * troll_q
  troll_log_pred <- log(troll_pred)
  troll_sig <- sqrt(troll_sd2)
  troll_log_obs <- log(troll_obs)
  troll_res <- (troll_obs - troll_pred) / troll_sig
  #if (troll_switch > 0) lp <- log(troll_sig) + 0.5 * troll_res^2
  if (troll_switch > 0) lp <- -dnorm(troll_log_obs, troll_log_pred, troll_sig, log = TRUE) 
  return(list(pred = troll_pred, resid = troll_res, lp = lp))
}

#' Genetic Tagging Likelihood
#'
#' Computes binomial likelihood of recapture events from genetic tagging data.
#'
#' @param gt_switch Integer switch to activate likelihood.
#' @param gt_obs Matrix of GT data [year_rel, age_rel, year_rec, nrel, nscan, nrec].
#' @param number_ysa 3D array [year, season, age] of numbers-at-age.
#' @return Negative log-likelihood vector per observation.
#' @export
#'
get_GT_like <- function(gt_switch, gt_obs, number_ysa) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_gt <- nrow(gt_obs)
  lp <- numeric(n_gt)
  nrec <- gt_obs[, 6] 
  # nrec <- OBS(nrec) # flag it for simulation
  gt_q <- 1  # hardcoded in C++; could be data/parameter if needed
  for (i in seq_len(n_gt)) {
    yrel <- gt_obs[i, 1] + 1  # release year index (1-based)
    arel <- gt_obs[i, 2] + 1  # release age index (1-based)
    # yrec <- gt_obs[i, 3] # not used here
    nrel <- gt_obs[i, 4]
    nscan <- gt_obs[i, 5]
    pgt <- nrel / (gt_q * number_ysa[yrel, 1, arel]) # Expected probability of recapture. Uses season 1 numbers-at-age.
    # Avoid log(0) and log(negative) by bounding pgt
    # pgt <- min(max(pgt, 1e-10), 1 - 1e-10)
    # Binomial log-likelihood
    # lp[i] <- nrec * log(pgt) + (nscan - nrec) * log(1 - pgt)
    if (gt_switch > 0) lp[i] <- -dbinom(nrec[i], nscan, pgt, log = TRUE)
  }
  return(lp)
}

#' Calculate recruitment
#'
#' Computes recruitment based on Beverton-Holt with depensation and log-normal deviations.
#'
#' @param y Year index.
#' @param sbio Spawning biomass.
#' @param B0 Unfished biomass.
#' @param alpha, beta Beverton-Holt parameters.
#' @param sigma_r Lognormal SD of recruitment deviations.
#' @param rdev_y Recruitment deviations.
#' @param sr_dep Depensation parameter (default 1e-10).
#' @return Recruitment value (numeric).
#' @export
get_recruitment <- function(y, sbio, B0, alpha, beta, sigma_r, rdev_y, sr_dep = 1e-10) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_year <- length(rdev_y)
  rec <- (alpha * sbio) / (beta + sbio) * (1 - exp(log(0.5) * sbio / (sr_dep * B0))) * exp(rdev_y[y] - 0.5 * sigma_r^2)
  # rec <- (alpha * sbio) / (beta + sbio) * (1 - exp(log(0.5) * sbio / (sr_dep * B0))) * exp(rdev_y[y])
  return(rec)
}


#' Harvest rate calculation
#'
#' Computes age-specific harvest rates and F by fishery.
#'
#' @param y, s Year and season index.
#' @param first_yr, first_yr_catch Model and catch start years.
#' @param removal_switch_f Vector of slice-switch flags.
#' @param catch_obs_ysf Observed catch by year-season-fishery.
#' @param number_ysa 3D array of numbers-at-age.
#' @param sel_fya Selectivity array.
#' @param weight_fya Weight-at-age.
#' @param af_sliced_ysfa Sliced numbers-at-age.
#' @return List with age-specific harvest rates, fishery F, and penalty term.
#' @export
get_harvest_rate <- function(y, s, first_yr, first_yr_catch, removal_switch_f,
                             catch_obs_ysf, number_ysa, sel_fya, weight_fya,
                             af_sliced_ysfa) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_fishery <- dim(weight_fya)[1]
  n_year <- dim(sel_fya)[2]
  n_age <- dim(sel_fya)[3]
  yy <- y - (first_yr_catch - first_yr)
  F_f <- numeric(n_fishery)
  h_rate_fa <- array(0, dim = c(n_fishery, n_age))
  for (f in seq_len(n_fishery)) {
    if (catch_obs_ysf[yy, s, f] > 0) {
      if (removal_switch_f[f] == 0) {
        Nsum <- sum(number_ysa[y, s,] * sel_fya[f, y,] * weight_fya[f, y,]) + 1e-6
        F_f[f] <- catch_obs_ysf[yy, s, f] / Nsum
        h_rate_fa[f,] <- F_f[f] * sel_fya[f,y,]
      } else if (removal_switch_f[f] == 1) {
        Nsum <- sum(af_sliced_ysfa[y, s, f,] * weight_fya[f, y,]) + 1e-6
        h_rate_fa[f,] <- (catch_obs_ysf[yy, s, f] * af_sliced_ysfa[y, s, f,]) / (number_ysa[y, s,] * Nsum)
      }
    }
  }
  h_rate_a <- colSums(h_rate_fa[1:6,])
  tmp <- posfun(x = 1 - sum(F_f), eps = 0.001)
  return(list(h_rate_a = h_rate_a, h_rate_fa = h_rate_fa, penalty = tmp$penalty))
}

#' Estimate Temporal Autocorrelation in Recruitment Deviations
#'
#' Calculates the autocorrelation coefficient (rho) for recruitment deviations.
#'
#' @param first_yr First model year.
#' @param last_yr Last model year.
#' @param rdev Vector of recruitment deviations.
#' @return Estimated autocorrelation (tau_ac2).
#' @export
#'
get_rho <- function(first_yr, last_yr, rdev) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  i1 <- 1965 - first_yr - 1
  i2 <- last_yr - first_yr - 5 - 1
  n  <- i2 - i1
  t1 <- rdev[(i1 + 1):(i1 + n)]
  t2 <- rdev[(i1 + 2):(i1 + n + 1)]
  t1m <- mean(t1)
  t2m <- mean(t2)
  t1_t1m <- t1 - t1m
  t2_t2m <- t2 - t2m
  tau_ac2 <- mean((t1 + 1 - t1m) * t2_t2m) / (sqrt(mean(t1_t1m^2)) * sqrt(mean(t2_t2m^2)))
  return(tau_ac2)
}

#' Initial Numbers and Beverton-Holt Parameters
#'
#' Computes the initial equilibrium numbers-at-age, unfished recruitment (R0), and Beverton-Holt stock-recruitment parameters.
#'
#' @param B0 Unfished spawning biomass.
#' @param steep Beverton-Holt steepness.
#' @param M_a Vector of natural mortality-at-age.
#' @param phi_ya Matrix [year, age] of spawning output-at-age.
#'
#' @return A list containing:
#' \describe{
#'   \item{Ninit}{Initial numbers-at-age (vector).}
#'   \item{R0}{Unfished recruitment (scalar).}
#'   \item{alpha}{BH alpha parameter.}
#'   \item{beta}{BH beta parameter.}
#' }
#' @export
#'
get_initial_numbers <- function(B0, steep, M_a, phi_ya) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_age <- length(M_a)
  rel_N <- numeric(n_age)
  rel_N[1] <- 1
  for (ia in 2:n_age) rel_N[ia] <- rel_N[ia - 1] * exp(-M_a[ia - 1])
  rel_N[n_age] <- rel_N[n_age] / (1 - exp(-M_a[n_age]))
  R0 <- B0 / sum(phi_ya[1,] * rel_N)
  alpha <- (4 * steep * R0) / (5 * steep - 1)
  beta  <- (B0 * (1 - steep)) / (5 * steep - 1)
  return(list(Ninit = R0 * rel_N, R0 = R0, alpha = alpha, beta = beta))
}

get_selectivity2 <- function(n_age, max_age, first_yr, first_yr_catch, 
                             sel_min_age_f, sel_max_age_f, sel_end_f, sel_change_year_fy,
                             par_log_sel) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_fishery <- nrow(sel_change_year_fy)
  n_year <- ncol(sel_change_year_fy)
  ymin <- first_yr_catch - first_yr + 1
  sel_fya <- array(0, dim = c(n_fishery, n_year, n_age))
  for (f in seq_len(n_fishery)) {
    amin <- sel_min_age_f[f] + 1
    amax <- sel_max_age_f[f] + 1
    ipar <- 1
    for (y in ymin:n_year) {
      if (sel_change_year_fy[f, y] != 0) {
        sel_tmp <- exp(par_log_sel[[f]][ipar,])
        ipar <- ipar + 1
        sel_fya[f, y, amin:amax] <- sel_tmp / mean(sel_tmp)
        if (as.logical(sel_end_f[f]) && amax < max_age) {
          for (a in (amax + 1):n_age) {
            sel_fya[f, y, a] <- sel_fya[f, y, amax]
          }
        }
      } else {
        sel_fya[f, y, ] <- sel_fya[f, y - 1, ]
      }
    }
  }
  return(sel_fya)
}

#' Compute Selectivity-at-Age Across Years
#'
#' Constructs the selectivity-at-age array by fishery and year, incorporating initial values and changes.
#'
#' @param n_age Total number of ages.
#' @param max_age Maximum model age.
#' @param first_yr Model start year.
#' @param first_yr_catch Vector of first catch years per fishery.
#' @param sel_end_f Logical vector indicating whether to extend the final selectivity across remaining ages.
#' @param sel_change_year_fy Matrix [fishery, year] indicating change years.
#' @param par_sels_init_i Vector of initial selectivity log-values.
#' @param par_sels_change_i Vector of changes in selectivity (log-space).
#' @return 3D array [fishery, year, age] of selectivity values.
#' @export
get_selectivity <- function(n_age, max_age, first_yr, first_yr_catch, 
                            sel_min_age_f, sel_max_age_f, sel_end_f, sel_change_year_fy,
                            par_sels_init_i, par_sels_change_i) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_fishery <- nrow(sel_change_year_fy)
  n_year <- ncol(sel_change_year_fy)
  sel_fya <- array(0, dim = c(n_fishery, n_year, n_age))
  ymin <- first_yr_catch - first_yr + 1
  ipar <- 1
  jpar <- 1
  
  for (f in seq_len(n_fishery)) {
    amin <- sel_min_age_f[f] + 1
    amax <- sel_max_age_f[f] + 1
    sel_tmp <- numeric(amax - amin + 1)
    # Initial selectivity for the first year with catch
    for (a in amin:amax) {
      sel_tmp[a - amin + 1] <- exp(par_sels_init_i[ipar])
      ipar <- ipar + 1
    }
    mean_tmp <- mean(sel_tmp)
    for (a in amin:amax) {
      sel_fya[f, ymin, a] <- sel_tmp[a - amin + 1] / mean_tmp
    }
    if (as.logical(sel_end_f[f]) && amax < max_age) {
      for (a in (amax + 1):n_age) {
        sel_fya[f, ymin, a] <- sel_fya[f, ymin, amax]
      }
    }
    # Selectivity in subsequent years
    for (y in (ymin + 1):n_year) {
      if (sel_change_year_fy[f, y] != 0) {
        for (a in amin:amax) {
          sel_tmp[a - amin + 1] <- sel_fya[f, y - 1, a] * exp(par_sels_change_i[jpar])
          jpar <- jpar + 1
        }
        mean_tmp <- mean(sel_tmp)
        for (a in amin:amax) {
          sel_fya[f, y, a] <- sel_tmp[a - amin + 1] / mean_tmp
        }
        if (as.logical(sel_end_f[f]) && amax < max_age) {
          for (a in (amax + 1):n_age) {
            sel_fya[f, y, a] <- sel_fya[f, y, amax]
          }
        }
      } else {
        sel_fya[f, y, ] <- sel_fya[f, y - 1, ]
      }
    }
  }
  return(sel_fya)
}


#' Spawning output-at-age
#'
#' Calculates spawning output per recruit at age and year using length-based maturity and fecundity scaling.
#'
#' @param log_psi Scalar fecundity scaling exponent.
#' @param length_m50 Length at 50% maturity.
#' @param length_m95 Length at 95% maturity.
#' @param length_mu_ysa 3D array [year, season, age] of expected length-at-age.
#' @param length_sd_a Vector of length standard deviations at age.
#' @param dl_yal 3D array [year, age, length_bin] of length distribution weights.
#' @return Matrix [year + 1, age] of normalized spawning output-at-age.
#' @export
#'
get_phi <- function(log_psi, length_m50, length_m95, length_mu_ysa, length_sd_a, dl_yal) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_age <- length(length_sd_a)
  n_year <- nrow(length_mu_ysa)
  n_bins <- 15
  phi_ya <- matrix(0, n_year + 1, n_age)
  psi <- exp(log_psi)
  for (iy in seq_len(n_year)) {
    phi_a <- numeric(n_age)
    for (ia in seq_len(n_age)) {
      lref <- length_mu_ysa[iy, 2, ia] # season 2 for phi
      sdref <- length_sd_a[ia]
      llq <- lref - 1.98 * sdref
      if (llq < 0) llq <- 0
      luq <- lref + 1.98 * sdref
      ldel <- (luq - llq) / (n_bins - 1)
      for (il in seq_len(n_bins)) {
        ltmp <- llq + (il - 1) * ldel
        matl <- 1 / (1 + 19^((length_m50 - ltmp) / (length_m95 - length_m50)))
        phil <- ltmp^(3 * psi) * matl
        phi_a[ia] <- phi_a[ia] + dl_yal[iy, ia, il] * phil
      }
    }
    #phi_ya[iy,] <- phi_a / max(phi_a)
    phi_ya[iy,] <- phi_a / phi_a[n_age] # don't use max/monotonically inc. fn.
  }
  phi_ya[n_year + 1,] <- phi_ya[n_year,]
  return(phi_ya)
}

#' Natural mortality-at-age
#'
#' Constructs a vector of M-at-age values using a declining early-age curve and late-age increase.
#'
#' @param min_age, max_age Minimum and maximum model age.
#' @param age_increase_M Age at which M begins to increase again.
#' @param m0 M at age-1 and 2.
#' @param m4 M at age-4 (controls slope of decline).
#' @param m10 M at age-10 (base for flat zone).
#' @param m30 M at age-30 (terminal age M).
#' @return Vector of M-at-age values.
#' @export
#'
get_M <- function(min_age, max_age, age_increase_M, m0, m4, m10, m30) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_age <- max_age - min_age + 1
  n_increase_M <- age_increase_M - min_age + 1
  M_slope <- log((m0 - m4) / (m0 - m10)) / log(0.33333)
  M_inc <- (m30 - m10) / (max_age - age_increase_M)
  M_age <- numeric(n_age)
  M_age[1:2] <- m0
  for (i in 3:10) {
    M_age[i] <- m0 - (m0 - m10) * ((i - 2) / 9)^M_slope
  }
  if (n_increase_M >= 10) {
    M_age[11:n_increase_M] <- m10
  }
  if (n_age > n_increase_M) {
    for (i in (n_increase_M + 1):(n_age - 1)) {
      M_age[i] <- M_age[i - 1] + M_inc
    }
  }
  M_age[n_age] <- m30
  return(M_age)
}


#' Transformation to ensure correlation parameters are between -1 and 1
#' @param x A numeric value.
#' @return A transformed value between -1 and 1.
rho_trans <- function(x) {
  2 / (1 + exp(-2 * x)) - 1
}

#' Positive Constraint Penalty Function
#' https://github.com/kaskr/adcomp/issues/7
#'
#' Ensures a value remains above a small threshold using a smooth approximation and penalty.
#'
#' @param x Numeric value to constrain.
#' @param eps Minimum allowable value (default 0.001).
#'
#' @return A list with:
#' \describe{
#'   \item{new}{Transformed value.}
#'   \item{penalty}{Penalty applied if \code{x < eps}.}
#' }
#' @importFrom RTMB ADoverload logspace_add
#' @export
#'
posfun <- function(x, eps = 0.001) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  pen <- eps * (1 / (1 - (x - eps) / eps + 
                       (x - eps)^2 / eps^2 - (x - eps)^3 / eps^3 + 
                       (x - eps)^4 / eps^4 - (x - eps)^5 / eps^5))
  out <- list()
  out$new <- eps * logspace_add(x / eps, 0)
  out$penalty <- pen
  return(out)
}



plot_cpue_lf <- function(data, object, posterior = NULL, probs = c(0.025, 0.975),
                    years = NULL, ...) {
  
  specs <- data.frame(Year = data$cpue_year + data$first_yr - 1, N = data$cpue_n) %>%
    mutate(id = 1:n())
  
  obs <- cbind(specs, data$cpue_lf) %>%
    data.frame() %>%
    pivot_longer(cols = starts_with("X"), names_to = "Length", values_to = "obs") %>%
    mutate(Length = parse_number(Length))
  
  pred <- cbind(specs, object$report()$cpue_lf_pred) %>%
    data.frame() %>%
    pivot_longer(cols = starts_with("X"), names_to = "Length", values_to = "pred") %>%
    mutate(Length = parse_number(Length))
  
  df <- full_join(obs, pred, by = join_by("Year", "N", "Length", "id")) %>% 
    mutate(Length = seq(87.5, 184, 4)[Length])

  if (!is.null(years)) df <- df %>% filter(Year %in% years)
  
  dfN <- df %>% 
    select(Year, N) %>%
    distinct() %>%
    mutate(N = paste0("N=", N))
  
  p <- ggplot(data = df, aes(x = .data$Length, y = .data$obs)) +
    geom_label(data = dfN, aes(x = -Inf, y = Inf, label = N), hjust = 0, vjust = 1, label.r = unit(0, "lines")) +
    geom_point(colour = "red") +
    geom_line(aes(y = .data$pred), linetype = "dashed") +
    labs(x = "Length (cm)", y = "Proportion") +
    facet_wrap(Year ~ .) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  
  # if (!is.null(posterior)) {
  #   df0 <- get_posterior(object = object, posterior = posterior, pars = "lf_pred") %>%
  #     mutate(Length = rep(1:25, each = 85)[id]) %>%
  #     mutate(id = rep(1:85, 31)[id]) %>%
  #     left_join(specs, by = join_by("id")) %>%
  #     select(-id, -output) %>%
  #     rename(pred = value) %>%
  #     filter(Fishery == fishery)
  #   
  #   df1 <- df0 %>% pivot_wider(names_from = Age, values_from = pred)
  #   prob <- as.matrix(df1 %>% select(`0`:`30`))
  #   
  #   df_ppred <- t(mapply(rmultinom, n = 1, size = df1$N, prob = split(x = prob, f = c(row(prob)))))
  #   df_ppred <- df_ppred / rowSums(df_ppred)
  #   
  #   dfpp <- cbind(df1 %>% select(chain, iter, Year, Fishery, N, min, max), df_ppred) %>%
  #     pivot_longer(cols = !chain:max, names_to = "Age", values_to = "ppred") %>%
  #     mutate(Age = as.numeric(Age) - 1)
  #   
  #   df_mcmc <- full_join(df0, dfpp, by = join_by("chain", "iter", "Age", "Year", "Fishery", "N", "min", "max")) %>%
  #     filter(Age >= min, Age <= max)
  #   
  #   p <- p +
  #     stat_summary(data = df_mcmc, geom = "ribbon", alpha = 0.5,
  #                  aes(y = ppred),
  #                  fun.min = function(x) quantile(x, probs = probs[1]),
  #                  fun.max = function(x) quantile(x, probs = probs[2])) +
  #     stat_summary(data = df_mcmc, geom = "ribbon", alpha = 0.5,
  #                  aes(y = pred),
  #                  fun.min = function(x) quantile(x, probs = probs[1]),
  #                  fun.max = function(x) quantile(x, probs = probs[2])) +
  #     stat_summary(data = df_mcmc, aes(y = pred), geom = "line", fun = median)
  # }
  
  return(p)
}
