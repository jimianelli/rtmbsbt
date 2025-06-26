#' The sbt model
#' 
#' Obtain the negative log-likelihood (NLL) value from the sbt model.
#' 
#' @param parameters a \code{list} of parameter values.
#' @param data a \code{list} of inputs.
#' @return the negative log-likelihood (NLL) value.
#' @importFrom RTMB getAll REPORT ADREPORT ADoverload plogis dnorm dlnorm dgamma dmultinom
#' @importFrom stats dcauchy
#' @export
#' 
sbt_model <- function(parameters, data) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  getAll(data, parameters, warn = FALSE)
  # aerial_log_obs <- OBS(aerial_log_obs)
  # cpue_log_obs <- OBS(cpue_log_obs)
  # troll_log_obs <- OBS(troll_log_obs)
  
  # Transformations
  
  par_B0 <- exp(par_log_B0)
  par_psi <- exp(par_log_psi)
  par_hsp_q <- exp(par_log_hsp_q)
  par_m0 <- exp(par_log_m0)
  par_m4 <- exp(par_log_m4)
  par_m10 <- exp(par_log_m10)
  par_m30 <- exp(par_log_m30)
  par_h <- exp(par_log_h)
  par_sigma_r <- exp(par_log_sigma_r)
  par_aerial_tau <- exp(par_log_aerial_tau)
  par_troll_tau <- exp(par_log_troll_tau)
  cpue_sigma <- exp(par_log_cpue_sigma)
  cpue_omega <- exp(par_log_cpue_omega)
  par_tag_H_factor <- exp(par_log_tag_H_factor)
  
  # State Variables
  
  M_a <- get_M(min_age, max_age, age_increase_M, par_m0, par_m4, par_m10, par_m30)
  S_a <- exp(-0.5 * M_a)
  phi_ya <- get_phi(par_psi, length_m50, length_m95, length_mu_ysa, length_sd_a, dl_yal)
  # for (i in 2:93) phi_ya[i,] <- phi_ya[1,]
  # sel_fya <- get_selectivity(n_age, max_age, first_yr, first_yr_catch,
  #                            sel_min_age_f, sel_max_age_f, sel_end_f,
  #                            sel_change_year_fy, par_sels_init_i, par_sels_change_i)
  par_log_sel_fya <- list(par_log_sel_1, par_log_sel_2, par_log_sel_3, par_log_sel_4, par_log_sel_5, par_log_sel_6)
  sel_fya <- get_selectivity2(n_age, max_age, first_yr, first_yr_catch,
                              sel_min_age_f, sel_max_age_f, sel_end_f,
                              sel_change_year_fy, par_log_sel_fya)
  lp_sel <- numeric(n_fishery)
  for (f in seq_len(n_fishery)) {
    f1 <- function(x) dautoreg(x, phi = exp(par_log_sel_phi[f, 1]), log = TRUE, scale = exp(par_log_sel_scale[f, 1])) # year
    f2 <- function(x) dautoreg(x, phi = exp(par_log_sel_phi[f, 2]), log = TRUE, scale = exp(par_log_sel_scale[f, 2])) # age
    lp_sel[f] <- -dseparable(f1, f2)(par_log_sel_fya[[f]])
  }
  
  # Initial conditions
  
  init <- get_initial_numbers(par_B0, par_h, M_a, phi_ya)
  R0 <- init$R0
  alpha <- init$alpha
  beta <- init$beta
  number_ysa <- array(0, dim = c(n_year + 1, n_season, n_age))
  number_ysa[1,1,] <- init$Ninit
  spawning_biomass_y <- numeric(n_year + 1)
  spawning_biomass_y[1] <- par_B0 # sum(number_ysa[1,1,] * phi_ya[1,])
  tau_ac2 <- get_rho(first_yr, last_yr, par_rdev_y)
  n_year2 <- n_year - 2
  rdev_y <- par_rdev_y
  for (y in n_year2:n_year) rdev_y[y] <- tau_ac2 * rdev_y[y - 1] + par_rdev_y[y]
  # rdev_y[] <- 0
  recruitment_y <- numeric(n_year)
  recruitment_y[1] <- R0
  
  hrate_ysa  <- array(0, dim = c(n_year + 1, n_season, n_age))
  hrate_ysfa  <- array(0, dim = c(n_year + 1, n_season, n_fishery, n_age))
  # F_ysf  <- array(0, dim = c(n_year + 1, n_season, n_fishery))
  catch_pred_fya <- array(0, dim = c(n_fishery, n_year + 1, n_age))
  catch_pred_ysf <- array(0, dim = c(n_year + 1, n_season, n_fishery))
  fy <- first_yr_catch - first_yr + 1
  n_age1 <- n_age - 1
  lp_penalty <- 0
  
  # Main population loop
  
  for (y in seq_len(n_year)) {
    # Season 1
    if (y >= fy) {
      hr <- get_harvest_rate(y, 1, first_yr, first_yr_catch, removal_switch_f, catch_obs_ysf, number_ysa, sel_fya, weight_fya, af_sliced_ysfa)
      hrate_ysa[y,1,] <- hr$h_rate_a
      hrate_ysfa[y,1,,] <- hr$h_rate_fa
      # F_ysf[y,1,] <- hr$F_f
      lp_penalty <- lp_penalty + hr$penalty
    }
    number_ysa[y,2,] <- number_ysa[y,1,] * (1 - hrate_ysa[y,1,]) * S_a
    # Season 2
    if (y >= fy) {
      hr <- get_harvest_rate(y, 2, first_yr, first_yr_catch, removal_switch_f, catch_obs_ysf, number_ysa, sel_fya, weight_fya, af_sliced_ysfa)
      hrate_ysa[y,2,] <- hr$h_rate_a
      hrate_ysfa[y,2,,] <- hr$h_rate_fa
      # F_ysf[y,2,] <- hr$F_f
      lp_penalty <- lp_penalty + hr$penalty
    }
    number_ysa[y + 1, 1, 2:n_age] <- number_ysa[y, 2, 1:n_age1] * (1 - hrate_ysa[y, 2, 1:n_age1]) * S_a[1:n_age1]
    number_ysa[y + 1, 1, n_age] <- number_ysa[y + 1, 1, n_age] + (number_ysa[y, 2, n_age] * (1 - hrate_ysa[y, 2, n_age]) * S_a[n_age])
    spawning_biomass_y[y + 1] <- sum(number_ysa[y + 1, 1,] * phi_ya[y + 1,])
    recruitment_y[y + 1] <- get_recruitment(y = y, sbio = spawning_biomass_y[y + 1], B0 = par_B0, alpha, beta, sigma_r = par_sigma_r, rdev_y)
    number_ysa[y + 1, 1, 1] <- recruitment_y[y + 1]
    for (f in seq_len(n_fishery)) {
      for (s in seq_len(n_season)) {
        catch_pred_fya[f, y,] <- catch_pred_fya[f, y,] + hrate_ysfa[y, s, f,] * number_ysa[y, s,]
        # catch_pred_fya[f, y,] <- catch_pred_fya[f, y,] + F_ysf[y, s, f] * sel_fya[f, y,] * number_ysa[y, s,]
        # for (a in seq_len(n_age)) {
        #   catch_pred_ysf[y, s, f] <- catch_pred_ysf[y, s, f] + F_ysf[y, s, f] * sel_fya[f, y, a] * number_ysa[y, s, a] * weight_fya[f, y, a]
      }
    }
  }

  # Priors
  
  lp_rec <- get_recruitment_prior(par_rdev_y, par_sigma_r, tau_ac2)
  lp_m10 <- 0
  lp_h <- 0
  lp_cpue_omega <- 0
  
  # Likelihoods
  
  x <- get_length_like(removal_switch_f, lf_year, lf_season, lf_fishery, lf_minbin, lf_obs, lf_n, catch_pred_fya, alk_ysal)
  lp_lf <- x$lp
  lf_pred <- x$pred
  x <- get_age_like(removal_switch_f, af_year, af_fishery, af_min_age, af_max_age, af_obs, af_n, catch_pred_fya)
  lp_af <- x$lp
  af_pred <- x$pred
  x <- get_cpue_like(cpue_switch, cpue_a1, cpue_a2, cpue_years, cpue_obs, cpue_adjust, cpue_sigma, cpue_omega, par_log_cpue_q, par_cpue_creep, number_ysa, sel_fya)
  lp_cpue <- x$lp
  cpue_pred <- x$pred
  cpue_resid <- x$resid
  x <- get_aerial_survey_like(aerial_switch, aerial_years, aerial_obs, aerial_cv, aerial_cov, par_aerial_tau, par_log_aerial_sel, number_ysa, weight_fya)
  lp_aerial <- x$lp
  lp_aerial_tau <- x$lp_aerial_tau
  aerial_pred <- x$pred
  aerial_resid <- x$resid
  x <- get_troll_like(troll_switch, troll_years, troll_obs, troll_sd, par_troll_tau, number_ysa)
  lp_troll <- x$lp
  troll_pred <- x$pred
  troll_resid <- x$resid
  x <- get_tag_like(tag_switch, min_K + 1, n_K, n_T, n_I, n_J, first_yr, M_a, hrate_ysa,
                    tag_release_cta, tag_recap_ctaa,
                    minI = tag_rel_min_age, maxI = tag_rel_max_age, maxJ = tag_recap_max_age,
                    shed1 = tag_shed_immediate, shed2 = tag_shed_continuous,
                    tag_rep_rates_ya, tag_H_factor = par_tag_H_factor, tag_var_factor)
  lp_tags <- x$lp
  tag_pred <- x$pred
  tag_resid <- x$resid
  lp_pop <- get_POP_like(pop_switch, pop_obs, phi_ya, paly, spawning_biomass_y)
  lp_hsp <- get_HSP_like(hsp_switch, hsp_obs, par_hsp_q, hsp_false_negative, number_ysa, phi_ya, M_a, spawning_biomass_y, hrate_ysa)
  lp_gt <- get_GT_like(gt_switch, gt_obs, number_ysa)
  
  nll <- sum(lp_sel) + lp_rec  + lp_m10 + lp_h + lp_cpue_omega +
    sum(lp_lf) + sum(lp_af) + sum(lp_cpue) + lp_aerial_tau + sum(lp_aerial) +
    sum(lp_troll) + sum(lp_tags) + sum(lp_pop) + sum(lp_hsp) + sum(lp_gt) + lp_penalty
  
  # Reporting
  
  REPORT(par_B0)
  REPORT(par_psi)
  REPORT(par_hsp_q)
  REPORT(par_m0)
  REPORT(par_m4)
  REPORT(par_m10)
  REPORT(par_m30)
  REPORT(par_h)
  REPORT(par_sigma_r)
  REPORT(cpue_sigma)
  REPORT(cpue_omega)
  
  REPORT(lp_sel)
  REPORT(lp_rec)
  REPORT(lp_aerial_tau)
  REPORT(lp_lf)
  REPORT(lp_af)
  REPORT(lp_cpue)
  REPORT(lp_aerial)
  REPORT(lp_troll)
  REPORT(lp_tags)
  REPORT(lp_pop)
  REPORT(lp_hsp)
  REPORT(lp_gt)
  REPORT(nll)
  
  REPORT(lf_pred)
  REPORT(af_pred)
  REPORT(cpue_pred)
  REPORT(cpue_resid)
  REPORT(aerial_pred)
  REPORT(aerial_resid)
  REPORT(troll_pred)
  REPORT(troll_resid)
  REPORT(tag_pred)
  REPORT(tag_resid)
  
  REPORT(M_a)
  REPORT(phi_ya)
  REPORT(spawning_biomass_y)
  REPORT(sel_fya)
  REPORT(number_ysa)
  REPORT(hrate_ysa)
  REPORT(catch_pred_ysf)
  REPORT(catch_pred_fya)
  
  REPORT(R0)
  REPORT(alpha)
  REPORT(beta)
  ADREPORT(par_sigma_r)
  REPORT(tau_ac2)
  REPORT(par_rdev_y)
  # REPORT(rec_dev_y)
  REPORT(rdev_y)
  REPORT(recruitment_y)
  
  return(nll)
}
