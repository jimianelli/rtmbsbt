

get_tag_like <- function(tag_switch, n_K, n_T, n_I, n_J, minK, M_a, hrate_ysa, 
                         par_hstar_i, tag_release_cta, tag_recap_ctaa, minI, maxI, maxJ, 
                         shed1, shed2, tag_rep_rates_ya,
                         tag_H_factor, tag_var_factor, tag_offset) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  lp <- 0
  hstar_s1_ya <- matrix(0, nrow = n_K, ncol = n_I)
  ipar <- 1
  for (k in seq_len(n_K)) {
    for (i in minI[k]:maxI[k]) {
      hstar_s1_ya[k, i] <- par_hstar_i[ipar]
      ipar <- ipar + 1
    }
  }
  
  # Arrays for results
  prR <- array(0, dim = c(n_K, n_T, n_I, n_J))
  tag_pred <- array(0, dim = c(n_K, n_T, n_I, n_J))
  tag_resid <- array(0, dim = c(n_K, n_T, n_I, n_J))
  
  # Survival and exploitation rates for tagged fish
  for (k in seq_len(n_K)) {
    for (j in minI[k]:maxJ[k]) {
      iy <- minK + k - 1 + j - 1
      hplus <- tag_H_factor * hrate_ysa[iy, 1, j]
      for (t in seq_len(n_T)) {
        # Survival/exploitation rates (season 1 and 2)
        S1 <- (1 - hplus) * (1 - hrate_ysa[iy, 2, j]) * exp(-M_a[j] - shed2[t])
        S2 <- (1 - hplus) * (1 - hrate_ysa[iy, 2, j]) * exp(-M_a[j] - 2 * shed2[t])
        f1 <- hplus + (1 - hplus) * exp(-0.5 * (M_a[j] + shed2[t])) * hrate_ysa[iy, 2, j]
        f2 <- hplus + (1 - hplus) * exp(-0.5 * (M_a[j] + 2 * shed2[t])) * hrate_ysa[iy, 2, j]
        # Star versions (using hstar_s1_ya)
        for (i in minI[k]:maxI[k]) {
          iy_star <- minK + k - 1 + i - 1
          S1star <- (1 - hstar_s1_ya[k, i]) * (1 - hrate_ysa[iy_star, 2, i]) * exp(-M_a[i] - shed2[t])
          S2star <- (1 - hstar_s1_ya[k, i]) * (1 - hrate_ysa[iy_star, 2, i]) * exp(-M_a[i] - 2 * shed2[t])
          f1star <- hstar_s1_ya[k, i] + (1 - hstar_s1_ya[k, i]) * exp(-0.5 * (M_a[i] + shed2[t])) * hrate_ysa[iy_star, 2, i]
          f2star <- hstar_s1_ya[k, i] + (1 - hstar_s1_ya[k, i]) * exp(-0.5 * (M_a[i] + 2 * shed2[t])) * hrate_ysa[iy_star, 2, i]
          for (j2 in minI[k]:maxJ[k]) {
            # Probabilities of recapture (prR)
            if (j2 < i) {
              prR[k, t, i, j2] <- 0
            } else if (j2 == i) {
              prR[k, t, i, j2] <- (2 * shed1[t] * f1star - shed1[t]^2 * f2star) * tag_rep_rates_ya[k + j2 - 2, j2]
            } else if (j2 == (i + 1)) {
              prR[k, t, i, j2] <- (2 * shed1[t] * S1star * f1 - shed1[t]^2 * S2star * f2) * tag_rep_rates_ya[k + j2 - 2, j2]
            } else if (j2 > (i + 1)) {
              prodS1 <- 1
              prodS2 <- 1
              for (s in (i + 1):(j2 - 1)) {
                prodS1 <- prodS1 * S1
                prodS2 <- prodS2 * S2
              }
              prR[k, t, i, j2] <- (2 * shed1[t] * S1star * prodS1 * f1 - shed1[t]^2 * S2star * prodS2 * f2) * tag_rep_rates_ya[k + j2 - 2, j2]
            }
          }
        }
      }
    }
  }
  
  # Predicted numbers of recaptures and residuals
  for (k in seq_len(n_K)) {
    for (t in seq_len(n_T)) {
      for (i in minI[k]:maxI[k]) {
        for (j in minI[k]:maxJ[k]) {
          tag_pred[k, t, i, j] <- tag_release_cta[k, t, i] * prR[k, t, i, j]
          tag_resid[k, t, i, j] <- (tag_recap_ctaa[k, t, i, j] - tag_pred[k, t, i, j]) /
            sqrt(1e-5 + tag_var_factor * tag_pred[k, t, i, j] * (1 - prR[k, t, i, j]))
        }
      }
    }
  }
  
  # Calculate likelihood
  if (tag_switch > 0) {
    for (k in seq_len(n_K)) {
      for (t in seq_len(n_T)) {
        for (i in minI[k]:maxI[k]) {
          totR <- 0
          totprR <- 0
          tag_od <- (tag_release_cta[k, t, i] - tag_var_factor) / (tag_var_factor - 1)
          if (tag_od < 0) tag_od <- 0.001
          lp <- lp + lgamma(tag_od) - lgamma(tag_release_cta[k, t, i] + tag_od)
          for (j in i:maxJ[k]) {
            lp <- lp + lgamma(tag_recap_ctaa[k, t, i, j] + tag_od * prR[k, t, i, j]) -
              lgamma(tag_od * prR[k, t, i, j])
            totR <- totR + tag_recap_ctaa[k, t, i, j]
            totprR <- totprR + prR[k, t, i, j]
          }
          notR <- tag_release_cta[k, t, i] - totR
          pr_notR <- 1 - totprR
          lp <- lp + lgamma(notR + tag_od * pr_notR) - lgamma(tag_od * pr_notR)
        }
      }
    }
    lp <- -lp + tag_offset
  } else {
    lp <- 0
  }
  
  # Returns negative log-likelihood as a vector of length 1
  return(lp)
}

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
  
  # Aerial selectivity (3 ages: 2, 3, 4)
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
  for (i in seq_len(n_aerial)) {
    y <- aerial_years[i]
    for (a in 1:3) { # Ages 2,3,4 = indices 3,4,5 (1-based)
      aerial_pred[i] <- aerial_pred[i] + number_ysa[y, 1, a + 2] * weight_fya[6, y, a + 2] * aerial_sel[a]
    }
  }
  aerial_resid <- log(aerial_obs) - log(aerial_pred) # Residuals (log obs - log pred)
  aerial_log_q <- sum(cov_inv %*% aerial_resid) / sum(cov_inv) # Estimate log_q
  aerial_resid <- aerial_resid - aerial_log_q
  aerial_pred <- aerial_pred * exp(aerial_log_q) # scale predictions
  lp_aerial_tau <- 0.5 * log(det(cov_matrix))
  if (aerial_switch > 0) {
    lp <- 0.5 * as.vector(aerial_resid %*% cov_inv %*% aerial_resid) # Mahalanobis term: 0.5 * x' Sigma^{-1} x
  } else {
    lp <- 0
  }
  return(list(pred = aerial_pred, resid = aerial_resid, lp_aerial_tau = lp_aerial_tau, lp = lp))
}

get_POP_like <- function(pop_switch, pop_obs, phi_ya, spawning_biomass_y) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_pops <- nrow(pop_obs)
  n_age <- ncol(phi_ya)
  lp <- numeric(n_pops)
  for (i in seq_len(n_pops)) {
    cc <- pop_obs[i, 1]
    ba <- pop_obs[i, 2] + 1
    nP <- pop_obs[i, 3]
    nC <- pop_obs[i, 4]
    pp <- (2 * phi_ya[cc, ba]) / spawning_biomass_y[cc] # parental probability
    # Avoid log(0)
    # pp <- min(max(pp, 1e-12), 1 - 1e-12) # NEED TO USE THE PENALTY CODE FROM CRA
    # Binomial log-likelihood
    # if (pop_switch > 0 && pp > 0) {
    if (nP > 0) {
      lp[i] <- -(nP * log(pp) + (nC - nP) * log(1 - pp))
    }# else {
    #   lp[i] <- -nC * log(1 - pp)
      # }
    # }
  }
  return(lp)
}

get_HSP_like <- function(hsp_switch, hsp_obs, hsp_q, hsp_false_negative, 
                         number_ysa, phi_ya, M_a, spawning_biomass_y, hrate_ysa) {

  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_hsp <- nrow(hsp_obs)
  n_year <- nrow(phi_ya)
  n_age <- ncol(phi_ya)
  lp <- numeric(n_hsp)
  
  for (i in seq_len(n_hsp)) {
    c1 <- hsp_obs[i, 1]
    c2 <- hsp_obs[i, 2]
    cdif <- hsp_obs[i, 3]
    nC <- hsp_obs[i, 4]
    nK <- hsp_obs[i, 5]
    cmin <- min(c1, c2)
    cmax <- max(c1, c2)
    xtmp <- 0
    for (a in seq_len(n_age)) {
      # Cumulative survival between c1 and c2 given reference adult age
      cumS <- 1
      if ((a + cdif) <= n_age) {
        for (ia in seq_len(cdif)) {
          age_idx <- a + ia - 1
          year_idx <- cmin + ia - 1
          cumS <- cumS * exp(-M_a[age_idx]) * (1 - hrate_ysa[year_idx, 1, age_idx]) * (1 - hrate_ysa[year_idx, 2, age_idx])
        }
      } else {
        for (ia in seq_len(cdif)) {
          age_idx <- a + ia - 1
          year_idx <- cmin + ia - 1
          idx <- if (age_idx <= n_age) age_idx else n_age
          cumS <- cumS * exp(-M_a[idx]) * (1 - hrate_ysa[year_idx, 1, idx]) * (1 - hrate_ysa[year_idx, 2, idx])
        }
      }
      # Effective increase in RO-at-age
      if ((a + cdif) <= n_age) {
        phi_val <- phi_ya[cmax, a + cdif]
      } else {
        phi_val <- phi_ya[cmax, n_age]
      }
      xtmp <- xtmp + number_ysa[cmin, 2, a] * phi_ya[cmin, a] / spawning_biomass_y[cmin] * cumS * phi_val
    }
    pp <- 4 * hsp_q * xtmp / spawning_biomass_y[cmax]
    phsp <- pp * hsp_false_negative
    
    # Likelihood calculation
    if (hsp_switch > 0 && phsp > 0) {
      if (nK > 0) {
        lp[i] <- -(nK * log(phsp) + (nC - nK) * log(1 - phsp))
      } else {
        lp[i] <- -nC * log(1 - phsp)
      }
    } else {
      lp[i] <- 0
    }
  }
  return(lp)
}


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

get_recruitment_prior <- function(rdev_y, sigma_r, tau_ac2) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_year <- length(rdev_y)
  r1 <- rdev_y[1:(n_year - 3)]
  r2 <- rdev_y[(n_year - 2):n_year]
  lp <- n_year * log(sigma_r)
  lp <- lp + 0.5 * sum(r1^2) / sigma_r^2
  lp <- lp + 0.5 * sum(r2^2) / (sigma_r^2 * (1 - tau_ac2^2))
  return(lp)
}


get_cpue_like <- function(cpue_switch, cpue_a1 = 5, cpue_a2 = 17, 
                          cpue_years, cpue_obs, cpue_adjust, cpue_sigma, cpue_omega, 
                          log_cpue_q, number_ysa, sel_fya) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  n_cpue <- length(cpue_obs)
  n_age <- dim(sel_fya)[3]
  cpue_log_pred <- numeric(n_cpue)
  for (i in seq_len(n_cpue)) {
    y <- cpue_years[i]
    # Get selectivity for the CPUE ages for this year (fishery 1, by convention)
    cpue_sel <- sel_fya[1, y, 5:n_age]
    cpue_n <- number_ysa[y, 2, 5:n_age]
    cpue_selm <- sel_fya[1, y, (cpue_a1 + 1):(cpue_a2 + 1)]
    tmpN <- sum(cpue_sel * cpue_n) / mean(cpue_selm)
    cpue_log_pred[i] <- log(cpue_adjust[i]) + cpue_omega * log(tmpN)
  }
  cpue_log_pred <- cpue_log_pred - log(mean(exp(cpue_log_pred))) + log_cpue_q
  cpue_pred <- exp(cpue_log_pred)
  cpue_resid <- (log(cpue_obs) - cpue_log_pred) / cpue_sigma
  if (cpue_switch > 0) {
    lp <- log(cpue_sigma) + 0.5 * cpue_resid^2
  } else {
    lp <- 0
  }
  return(list(pred = cpue_pred, resid = cpue_resid, lp = lp))
}

get_age_like <- function(af_year, af_season, af_fishery, af_minage, af_obs, af_n, catch_pred_fya) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  n_af <- nrow(af_obs)
  n_bins <- ncol(af_obs)
  n_age <- dim(catch_pred_fya)[3]
  lp <- numeric(n_af)
  age_pred <- matrix(0, n_af, n_bins)
  
  for (i in seq_len(n_af)) {
    y <- af_year[i]
    s <- af_season[i]
    f <- af_fishery[i]
    minage <- af_minage[i]
    pred_age <- numeric(n_bins)
    for (a in seq_len(n_bins)) {
      pred_age[a] <- catch_pred_fya[f, y, a + minage - 1]
    }
    total_pred <- sum(pred_age)
    if (total_pred > 0) {
      pred_age <- pred_age / total_pred
    } else {
      pred_age <- 1 / n_bins
    }
    age_pred[i, ] <- pred_age
    # Multinomial negative log-likelihood
    obs <- af_obs[i, ]
    # Add small value to avoid log(0)
    pred_age <- pmax(pred_age, 1e-10)
    lp[i] <- -sum(obs * log(pred_age))
    # Optionally, you could use full multinomial likelihood:
    # lp[i] <- - (lgamma(af_n[i] + 1) - sum(lgamma(obs + 1)) + sum(obs * log(pred_age)))
  }
  return(lp)
}


get_length_like <- function(lf_year, lf_season, lf_fishery, lf_minbin, lf_obs, 
                            lf_n, catch_pred_fya, alk_ysal) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  n_lf <- nrow(lf_obs)
  n_bins <- ncol(lf_obs)
  n_age <- dim(catch_pred_fya)[3]
  lp <- numeric(n_lf)
  lf_pred <- matrix(0, n_lf, n_bins)
  
  for (i in seq_len(n_lf)) {
    y <- lf_year[i]
    s <- lf_season[i]
    f <- lf_fishery[i]
    minbin <- lf_minbin[i]
    pred_len <- numeric(n_bins)
    # Predict length freq by summing predicted catch at age * ALK for each bin
    for (a in seq_len(n_age)) {
      # for (l in seq_len(n_bins)) {
        # pred_len[l] <- pred_len[l] + catch_pred_fya[f, y, a] * alk_ysal[y, s, a, l + minbin - 1]
        pred_len[l] <- pred_len[l] + catch_pred_fya[f, y, a] * alk_ysal[y, s, a,]
      # }
    }
    total_pred <- sum(pred_len)
    if (total_pred > 0) {
      pred_len <- pred_len / total_pred
    } else {
      pred_len[] <- 1 / n_bins
    }
    lf_pred[i, ] <- pred_len
    # Multinomial log-likelihood
    obs <- lf_obs[i, ]
    n <- lf_n[i]
    # Add small value to avoid log(0)
    pred_len <- pmax(pred_len, 1e-10)
    lp[i] <- -sum(obs * log(pred_len))
    # Optionally, add log multinomial coefficient if modeling full multinomial:
    # lp[i] <- - (lgamma(n + 1) - sum(lgamma(obs + 1)) + sum(obs * log(pred_len)))
  }
  return(lp)
}

get_troll_like <- function(troll_switch, troll_years, troll_obs, troll_sd, par_troll_tau, number_ysa) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  # THIS ALL LOOKS WRONG
  n_troll <- length(troll_obs)
  troll_pred <- troll_resid <- lp <- numeric(n_troll)
  for (i in seq_len(n_troll)) {
    y <- troll_years[i]
    pred <- log(number_ysa[y, 1, 2]) + par_troll_tau
    troll_pred[i] <- pred
    troll_resid[i] <- troll_obs[i] - pred
    lp[i] <- 0.5 * (troll_resid[i]^2 / (troll_sd[i]^2)) + log(troll_sd[i]) + 0.5 * log(2 * pi)
  }
  if (troll_switch > 0) {
    return(lp)
  } else {
    return(rep(0, n_troll))
  }
}

get_GT_like <- function(gt_switch, gt_obs, number_ysa) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_gt <- nrow(gt_obs)
  lp <- numeric(n_gt)
  gt_q <- 1  # hardcoded in C++; could be data/parameter if needed
  for (i in seq_len(n_gt)) {
    yrel <- gt_obs[i, 1]   # release year index (1-based)
    arel <- gt_obs[i, 2]   # release age index (1-based)
    # yrec <- gt_obs[i, 3] # not used here
    nrel <- gt_obs[i, 4]
    nscan <- gt_obs[i, 5]
    nrec <- gt_obs[i, 6]
    # Expected probability of recapture. Uses season 1 numbers-at-age.
    pgt <- nrel / (gt_q * number_ysa[yrel, 1, arel])
    # Avoid log(0) and log(negative) by bounding pgt
    # pgt <- min(max(pgt, 1e-10), 1 - 1e-10)
    # Binomial log-likelihood
    lp[i] <- nrec * log(pgt) + (nscan - nrec) * log(1 - pgt)
  }
  
  if (gt_switch > 0) {
    return(-lp)
  } else {
    return(rep(0, n_gt))
  }
}

get_recruitment <- function(y, sbio, B0, alpha, beta, sigma_r, rdev_y) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_year <- length(rdev_y)
  sr_dep <- 1e-10 # depensation parameter
  # rec <- (alpha * sbio) / (beta + sbio) * (1 - exp(log(0.5) * sbio / (sr_dep * B0))) * exp(rdev_y[y] - 0.5 * sigma_r^2)
  rec <- (alpha * sbio) / (beta + sbio) * (1 - exp(log(0.5) * sbio / (sr_dep * B0))) * exp(rdev_y[y] - 0.5 * sigma_r^2)
  # rec <- (alpha * sbio) / (beta + sbio) * (1 - exp(log(0.5) * sbio / (sr_dep * B0))) * exp(rdev_y[y])
  return(rec)
}


get_harvest_rate <- function(y, s, first_yr, first_yr_catch, catch_obs_ysf, number_ysa, sel_fya, weight_fya) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  n_fishery <- dim(sel_fya)[1]
  n_year <- dim(sel_fya)[2]
  n_age <- dim(sel_fya)[3]
  yy <- y - (first_yr_catch - first_yr)
  F_f <- numeric(n_fishery)
  h_rate_a <- numeric(n_age)
  
  for (f in seq_len(n_fishery)) {
    if (catch_obs_ysf[yy, s, f] > 0) {
      Nsum <- 1e-6 + sum(number_ysa[y, s,] * sel_fya[f, y,] * weight_fya[f, y,])
      F_f[f] <- catch_obs_ysf[yy, s, f] / Nsum
      h_rate_a <- h_rate_a + F_f[f] * sel_fya[f,y,]
    }
  }
  # kap <- 100
  # for (a in seq_len(n_age)) {
  #   if (h_rate_a[a] > 0.9) {
  #     h_rate_a[a] <- 0.9
  #   }
  # }
  return(list(h_rate_a = h_rate_a, F_f = F_f))
}

get_rho <- function(first_yr, last_yr, rdev) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  # rdev: vector of recruitment deviations (length at least last_yr - first_yr + 1)
  # Returns: tau_ac2 (autocorrelation parameter)
  # Model years: 1931-2022
  # Rec years: 1932-2023
  # n_years = n_recs: 92
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

get_initial_numbers <- function(B0, steep, M_a, phi_ya) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  # B0: equilibrium (unfished) spawning biomass (scalar)
  # steep: steepness parameter (scalar)
  # M_a: vector of natural mortality at age (length n_age)
  # phi_ya: list of length n_year + 1, each element is a vector of length n_age (use phi_ya[[1]] for initial year)
  # Returns:
  #   rel_N: vector of initial numbers-at-age (length n_age)
  #   R0: equilibrium recruitment (scalar)
  #   alpha, beta: stock-recruitment parameters (scalars)
  n_age <- length(M_a)
  rel_N <- numeric(n_age)
  rel_N[1] <- 1
  for (ia in 2:n_age) rel_N[ia] <- rel_N[ia - 1] * exp(-M_a[ia - 1])
  rel_N[n_age] <- rel_N[n_age] / (1 - exp(-M_a[n_age]))
  R0 <- B0 / sum(phi_ya[1,] * rel_N)
  alpha <- (4 * steep * R0) / (5 * steep - 1)
  beta  <- (B0 * (1 - steep)) / (5 * steep - 1)
  list(Ninit = R0 * rel_N, R0 = R0, alpha = alpha, beta = beta)
}

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


get_phi <- function(psi, length_m50, length_m95, length_mu_ysa, length_sd_a, dl_yal) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  # psi, length_m50, length_m95: scalars
  # length_mu_ysa: [year, season, age] array (use season 2, i.e., [,2,] in R, since C++ uses 1-based and season=1 for phi)
  # length_sd_a: vector of length n_age
  # dl_yal: [year, age, bin] array
  # Returns: list of phi vectors, one for each year (length n_year + 1), each vector of length n_age
  n_age <- length(length_sd_a)
  n_year <- nrow(length_mu_ysa)
  n_bins <- 15
  phi_ya <- matrix(0, n_year + 1, n_age)
  
  for (iy in seq_len(n_year)) { # R is 1-based
    phi_a <- numeric(n_age)
    for (ia in seq_len(n_age)) {
      lref <- length_mu_ysa[iy, 2, ia] # season 2 for phi (C++ uses (iy,1,ia); check your data if indexing differs)
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

get_M <- function(min_age, max_age, age_increase_M, m0, m4, m10, m30) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_age <- max_age - min_age + 1
  n_increase_M <- age_increase_M - min_age + 1
  M_slope <- log((m0 - m4) / (m0 - m10)) / log(0.33333)
  M_inc <- (m30 - m10) / (max_age - age_increase_M)
  M_age <- numeric(n_age)
  M_age[1] <- m0
  M_age[2] <- m0
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
