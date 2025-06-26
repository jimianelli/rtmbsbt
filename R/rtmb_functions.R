#' Helper to make closure
#' 
#' @param f a \code{vector} of midpoints.
#' @param d natural mortality at the reference size.
#' @return a \code{vector}.
#' @export
#' 
cmb <- function(f, d) function(p) f(p, d)


get_tag_like <- function(tag_switch, minK, n_K, n_T, n_I, n_J, 
                         first_yr, M_a, hrate_ysa, 
                         par_hstar_i, tag_release_cta, tag_recap_ctaa, 
                         minI, maxI, maxJ, 
                         shed1, shed2, tag_rep_rates_ya,
                         tag_H_factor, tag_var_factor, tag_offset) {
  
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  hstar_s1_ya <- matrix(0, nrow = n_K, ncol = n_I)
  ipar <- 1
  for (k in seq_len(n_K)) {
    for (i in minI[k]:maxI[k]) {
      hstar_s1_ya[k, i] <- par_hstar_i[ipar]
      ipar <- ipar + 1
    }
  }
  
  prR <- tag_pred <- tag_resid <- array(0, dim = c(n_K, n_T, n_I, n_J))
  S1 <- S2 <- f1 <- f2 <- array(0, dim = c(n_K, n_T, n_J))
  S1star <- S2star <- f1star <- f2star <- array(0, dim = c(n_K, n_T, n_I))
  
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
        S1star[k, t, i] <- (1 - hstar_s1_ya[k, i]) * (1 - hrate_ysa[iy, 2, i + 1]) * exp(-M_a[i + 1] - shed2[t])
        S2star[k, t, i] <- (1 - hstar_s1_ya[k, i]) * (1 - hrate_ysa[iy, 2, i + 1]) * exp(-M_a[i + 1] - 2 * shed2[t])
        f1star[k, t, i] <- hstar_s1_ya[k, i] + (1 - hstar_s1_ya[k, i]) * exp(-0.5 * (M_a[i + 1] + shed2[t])) * hrate_ysa[iy, 2, i + 1]
        f2star[k, t, i] <- hstar_s1_ya[k, i] + (1 - hstar_s1_ya[k, i]) * exp(-0.5 * (M_a[i + 1] + 2 * shed2[t])) * hrate_ysa[iy, 2, i + 1]
      }
    }
  }
  
  # Generate probabilities of recapture
  # tag_rep_rates_ya: years 1991-1997, ages 1-8
  for (k in seq_len(n_K)) {
    for (t in seq_len(n_T)) {
      for (i in minI[k]:maxI[k]) {
        for (j in minI[k]:maxJ[k]) {
          if (j == i) {
            prR[k, t, i, j] <- (2 * shed1[t] * f1star[k, t, j] - shed1[t]^2 * f2star[k, t, j]) * tag_rep_rates_ya[k + j - 2, j]
          }
          if (j == (i + 1)) {
            prR[k, t, i, j] <- (2 * shed1[t] * S1star[k, t, i] * f1[k, t, j] - shed1[t]^2 * S2star[k, t, i] * f2[k, t, j]) * tag_rep_rates_ya[k + j - 2, j]
          }
          if (j > (i + 1)) {
            prodS1 <- 1
            prodS2 <- 1
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
        for (j in minI[k]:maxJ[k]) {
          tag_pred[k, t, i, j] <- tag_release_cta[k, t, i] * prR[k, t, i, j]
          tag_resid[k, t, i, j] <- (tag_recap_ctaa[k, t, i, j] - tag_pred[k, t, i, j]) / 
            sqrt(1e-5 + tag_var_factor * tag_pred[k, t, i, j] * (1 - prR[k, t, i, j]))
        }
      }
    }
  }
  
  # 7) preliminary overdispersion offset
  tag_offset1 <- 0
  for (k in seq_len(n_K)) {
    for (t in seq_len(n_T)) {
      for (i in minI[k]:maxI[k]) {
        Nrel <- tag_release_cta[k,t,i]
        tag_od <- (Nrel - tag_var_factor) / (tag_var_factor - 1)
        if (tag_od <= 0) tag_od <- 1e-3
        tag_offset1 <- tag_offset1 + lgamma(tag_od) - lgamma(Nrel + tag_od)
        totR <- 0 
        totprR <- 0
        for (j in i:maxJ[k]) {
          prR1 <- 1e-6 + tag_recap_ctaa[k,t,i,j] / (1e-6 + Nrel)
          tag_offset1 <- tag_offset1 + lgamma(tag_recap_ctaa[k,t,i,j] + tag_od * prR1) - lgamma(tag_od * prR1)
          totR   <- totR   + tag_recap_ctaa[k,t,i,j]
          totprR <- totprR + prR1
        }
        notR    <- Nrel - totR
        pr_notR <- 1 - totprR
        tag_offset1 <- tag_offset1 + lgamma(notR + tag_od * pr_notR) - lgamma(tag_od * pr_notR)
      }
    }
  }
  
  # 8) log‐likelihood under beta‐binomial
  lp <- 0
  if (tag_switch > 0) {
    loglkhd_R <- 0
    for (k in seq_len(n_K)) {
      for (t in seq_len(n_T)) {
        for (i in minI[k]:maxI[k]) {
          Nrel <- tag_release_cta[k,t,i]
          tag_od <- (Nrel - tag_var_factor) / (tag_var_factor - 1)
          if (tag_od < 0) tag_od <- 1e-3
          loglkhd_R <- loglkhd_R + lgamma(tag_od) - lgamma(Nrel + tag_od)
          totR <- 0
          totprR <- 0
          for (j in i:maxJ[k]) {
            loglkhd_R <- loglkhd_R + lgamma(tag_recap_ctaa[k,t,i,j] + tag_od * prR[k,t,i,j]) - lgamma(tag_od * prR[k,t,i,j]) 
            totR <- totR + tag_recap_ctaa[k,t,i,j]
            totprR <- totprR + prR[k,t,i,j]
          }
          notR <- Nrel - totR
          pr_notR <- 1 - totprR
          loglkhd_R <- loglkhd_R + lgamma(notR + tag_od * pr_notR) - lgamma(tag_od * pr_notR)
        }
      }
    }
    lp <- -loglkhd_R + tag_offset1
  }
  
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


get_POP_like <- function(pop_switch, pop_obs, phi_ya, paly, spawning_biomass_y) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_pops <- nrow(pop_obs)
  n_age <- ncol(phi_ya)
  lp <- numeric(n_pops)
  for (i in seq_len(n_pops)) {

    # adult has a direct age

    if(pop_obs[i,4] == 0) {

      cc <- pop_obs[i, 1] 
      yy <- pop_obs[i, 2]
      aa <- pop_obs[i, 3] + 1 # estimated age of adult @ time of capture
      ba <- ifelse(aa-(yy-cc) < 1,1,aa-(yy-cc))
      ba <- min(ba,n_age)
      nP <- pop_obs[i, 5]
      nC <- pop_obs[i, 6]
      pp <- (2 * phi_ya[cc, ba]) / spawning_biomass_y[cc] # parental probability

    }

    # adult has observed length only

    if (pop_obs[i,4] == 1) {
      cc <- pop_obs[i, 1] 
      yy <- pop_obs[i, 2]
      ll <- pop_obs[i, 3] # observed length bin of adult @ time of capture
      nP <- pop_obs[i, 5]
      nC <- pop_obs[i, 6]
      amin <- yy - cc + 1 # anything younger than this can't be a parent
      arng <- amin:n_age
      ba <- arng - (yy - cc)
      pp <- (2 / spawning_biomass_y[cc]) * sum(paly[ll,arng,yy] * phi_ya[cc, ba]) 
    }
    lp[i] <- -(nP * log(pp) + (nC - nP) * log(1 - pp)) 
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
  # Calculate the relative age distribution in the adults for each year and age
  gamx_ya <- matrix(0, nrow = n_year, ncol = n_age)
  for (y in seq_len(n_year)) {
    # Use season 2 for numbers-at-age of adults
    gamx_ya[y,] <- number_ysa[y, 2,] * phi_ya[y,] / spawning_biomass_y[y]
  }
  for (i in seq_len(n_hsp)) {
    cmin <- hsp_obs[i, 1] + 1
    cmax <- hsp_obs[i, 2] + 1
    cdif <- hsp_obs[i, 3]
    nC <- hsp_obs[i, 4]
    nK <- hsp_obs[i, 5]
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
          idx <- if (age_idx < n_age) age_idx else n_age
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
    
    # Likelihood calculation
    # if (hsp_switch > 0 && phsp > 0) {
      # if (nK > 0) {
        lp[i] <- -(nK * log(phsp) + (nC - nK) * log(1 - phsp))
      # } else {
      #   lp[i] <- -nC * log(1 - phsp)
      # }
    # } else {
    #   lp[i] <- 0
    # }
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
    lp <- numeric(n_cpue)
  }
  return(list(pred = cpue_pred, resid = cpue_resid, lp = lp))
}

get_age_like <- function(af_year, af_fishery, af_min_age, af_max_age, af_obs, af_n, catch_pred_fya) {
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
    lp[i] <- af_n[i] * sum(obs * log(obs)) - af_n[i] * sum(obs * log(pred))
  }
  return(list(pred = age_pred, lp = lp))
}


get_length_like <- function(lf_year, lf_season, lf_fishery, lf_minbin, lf_obs, 
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
      # pred[1:(mbin - 1)] <- 0
      obs[mbin] <- sum(obs[1:mbin])
      # obs[1:(mbin - 1)] <- 0
      obs <- obs[mbin:n_bins]
      pred <- pred[mbin:n_bins]
    }
    lf_pred[i, mbin:n_bins] <- pred
    lp[i] <- -lf_n[i] * sum(obs * log(pred))
    obs <- obs + 1e-6
    lp[i] <- lp[i] + lf_n[i] * sum(obs * log(obs))
  }
  return(list(pred = lf_pred, lp = lp))
}

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
  troll_sig <- sqrt(troll_sd2)
  troll_res <- (troll_obs - troll_pred) / troll_sig
  if (troll_switch > 0) lp <- log(troll_sig) + 0.5 * troll_res^2
  return(list(pred = troll_pred, resid = troll_res, lp = lp))
}

get_GT_like <- function(gt_switch, gt_obs, number_ysa) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_gt <- nrow(gt_obs)
  lp <- numeric(n_gt)
  gt_q <- 1  # hardcoded in C++; could be data/parameter if needed
  for (i in seq_len(n_gt)) {
    yrel <- gt_obs[i, 1] + 1  # release year index (1-based)
    arel <- gt_obs[i, 2] + 1  # release age index (1-based)
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

# depensation parameter
get_recruitment <- function(y, sbio, B0, alpha, beta, sigma_r, rdev_y, sr_dep = 1e-10) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_year <- length(rdev_y)
  rec <- (alpha * sbio) / (beta + sbio) * (1 - exp(log(0.5) * sbio / (sr_dep * B0))) * exp(rdev_y[y] - 0.5 * sigma_r^2)
  # rec <- (alpha * sbio) / (beta + sbio) * (1 - exp(log(0.5) * sbio / (sr_dep * B0))) * exp(rdev_y[y])
  return(rec)
}


get_harvest_rate <- function(y, s, first_yr, first_yr_catch, slice_switch_f,
                             catch_obs_ysf, number_ysa, sel_fya, weight_fya,
                             sliced_ysfa) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  n_fishery <- dim(sel_fya)[1]
  n_year <- dim(sel_fya)[2]
  n_age <- dim(sel_fya)[3]
  yy <- y - (first_yr_catch - first_yr)
  F_f <- numeric(n_fishery)
  h_rate_fa <- array(0, dim = c(n_fishery, n_age))

  for (f in seq_len(n_fishery)) {
    if (catch_obs_ysf[yy, s, f] > 0) {
      if (slice_switch_f[f] == 1) {
        Nsum <- sum(sliced_ysfa[y,s,f,] * weight_fya[f,y,]) + 1e-6
        h_rate_fa[f,] <- (catch_obs_ysf[y,s,f] * sliced_ysfa[y,s,f,]) / (number_ysa[y,s,] * Nsum)
      } else {
        Nsum <- sum(number_ysa[y, s,] * sel_fya[f, y,] * weight_fya[f, y,]) + 1e-6
        F_f[f] <- catch_obs_ysf[yy, s, f] / Nsum
        h_rate_fa[f,] <- F_f[f] * sel_fya[f,y,]
      }
    }
  }
  h_rate_a <- colSums(h_rate_fa)
  tmp <- posfun(x = 1 - sum(F_f), eps = 0.001)
  return(list(h_rate_a = h_rate_a, F_f = F_f, penalty = tmp$penalty))
}

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
  n_age <- length(length_sd_a)
  n_year <- nrow(length_mu_ysa)
  n_bins <- 15
  phi_ya <- matrix(0, n_year + 1, n_age)
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

#' Inverse logit transformation with bounds
#' @param x A numeric value.
#' @param lb Lower bound (default is 0).
#' @param ub Upper bound (default is 1).
#' @return A transformed value between lb and ub.
invlogit_bnd <- function(x, lb = 0, ub = 1) {
  expx <- exp(x)
  (expx + lb) / (ub + expx)
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

#' Ensure population above 0Add commentMore actions
#' 
#' https://github.com/kaskr/adcomp/issues/7
#' 
#' @param x value to remain above eps
#' @param eps value to compare x to
#' @return a \code{vector} penalty
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
