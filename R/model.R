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
  # catch_obs <- OBS(catch_obs)
  # cpue_obs <- OBS(cpue_obs)
  # lf_obs <- OBS(lf_obs) # doesn't work for multinomial without rounding and all sorts
  
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
  par_hstar_i <- plogis(par_logit_hstar_i)
  par_tag_H_factor <- exp(par_log_tag_H_factor)

  # State Variables
  
  M_a <- get_M(min_age, max_age, age_increase_M, par_m0, par_m4, par_m10, par_m30)
  S_a <- exp(-0.5 * M_a)
  phi_ya <- get_phi(par_psi, length_m50, length_m95, length_mu_ysa, length_sd_a, dl_yal)
  
  sel_fya <- get_selectivity(n_age, max_age, first_yr, first_yr_catch, 
                             sel_min_age_f, sel_max_age_f, sel_end_f, 
                             sel_change_year_fy, par_sels_init_i, par_sels_change_i)
  
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
  recruitment_y <- numeric(n_year)
  recruitment_y[1] <- R0
  
  # Main population loop
  hrate_ysa  <- array(0, dim = c(n_year + 1, n_season, n_age))
  F_ysf  <- array(0, dim = c(n_year + 1, n_season, n_fishery))
  catch_pred_fya <- array(0, dim = c(n_fishery, n_year + 1, n_age))
  catch_pred_ysf <- array(0, dim = c(n_year + 1, n_season, n_fishery))
  fy <- first_yr_catch - first_yr + 1
  n_age1 <- n_age - 1
  lp_penalty <- 0
  
  slice_switch_f <- numeric(n_fishery)
  sliced_ysfa <- 1
  
  for (y in seq_len(n_year)) {
    # Season 1
    if (y >= fy) {
      hr <- get_harvest_rate(y, 1, first_yr, first_yr_catch, slice_switch_f, catch_obs_ysf, number_ysa, sel_fya, weight_fya, sliced_ysfa)
      hrate_ysa[y,1,] <- hr$h_rate_a
      F_ysf[y,1,] <- hr$F_f
      lp_penalty <- lp_penalty + hr$penalty
    }
    number_ysa[y,2,] <- number_ysa[y,1,] * (1 - hrate_ysa[y,1,]) * S_a
    # Season 2
    if (y >= fy) {
      hr <- get_harvest_rate(y, 2, first_yr, first_yr_catch, slice_switch_f, catch_obs_ysf, number_ysa, sel_fya, weight_fya, sliced_ysfa)
      hrate_ysa[y,2,] <- hr$h_rate_a
      F_ysf[y,2,] <- hr$F_f
      lp_penalty <- lp_penalty + hr$penalty
    }
    number_ysa[y + 1, 1, 2:n_age] <- number_ysa[y, 2, 1:n_age1] * (1 - hrate_ysa[y, 2, 1:n_age1]) * S_a[1:n_age1]
    number_ysa[y + 1, 1, n_age] <- number_ysa[y + 1, 1, n_age] + (number_ysa[y, 2, n_age] * (1 - hrate_ysa[y, 2, n_age]) * S_a[n_age])
    spawning_biomass_y[y + 1] <- sum(number_ysa[y + 1, 1, ] * phi_ya[y + 1,])
    recruitment_y[y + 1] <- get_recruitment(y = y, sbio = spawning_biomass_y[y + 1], B0 = par_B0, alpha, beta, sigma_r = par_sigma_r, rdev_y)
    number_ysa[y + 1, 1, 1] <- recruitment_y[y + 1]
    for (f in seq_len(n_fishery)) {
      for (s in seq_len(n_season)) {
        catch_pred_fya[f, y,] <- catch_pred_fya[f, y,] + F_ysf[y, s, f] * sel_fya[f, y,] * number_ysa[y, s,]
        # for (a in seq_len(n_age)) {
        #   catch_pred_ysf[y, s, f] <- catch_pred_ysf[y, s, f] + F_ysf[y, s, f] * sel_fya[f, y, a] * number_ysa[y, s, a] * weight_fya[f, y, a]
      }
    }
  }

  # Likelihoods and priors
  
  lp_sel <- get_sel_like(first_yr, first_yr_catch_f, sel_min_age_f, sel_max_age_f, 
                         sel_change_year_fy, sel_change_sd_fy, sel_smooth_sd_f, 
                         par_sels_init_i, par_sels_change_i, sel_fya)
  
  # length(par_sels_init_i)
  # length(par_sels_change_i)
  # 
  # f1 <- function(x) dautoreg(x, phi = 0.3, log = TRUE)
  # f2 <- function(x) dautoreg(x, phi = 0.1, log = TRUE)
  # x <- array(rnorm(100), dim = c(3, 2))
  # -dseparable(f1, f2)(x)
  
  lp_rec <- get_recruitment_prior(par_rdev_y, par_sigma_r, tau_ac2)
  lp_hstar <- 0.1 * sum((log(par_hstar_i) + 6)^2)
  lp_m10 <- 0
  lp_h <- 0
  lp_cpue_omega <- 0

  # Data likelihoods
  
  lp_lf <- get_length_like(lf_year, lf_season, lf_fishery, lf_minbin, lf_obs, 
                           lf_n, catch_pred_fya, alk_ysal)
  # lf_pred <- matrix(0, n_lf, 25)
  lp_af <- get_age_like(af_year, af_fishery, af_min_age, af_max_age, af_obs, af_n, catch_pred_fya)
  af_pred <- matrix(0, n_af, n_age)
  x <- get_cpue_like(cpue_switch, cpue_a1, cpue_a2, cpue_years, cpue_obs, cpue_adjust, cpue_sigma, cpue_omega, par_log_cpue_q, number_ysa, sel_fya)
  lp_cpue <- x$lp
  cpue_pred <- x$pred
  cpue_resid <- x$resid
  x <- get_aerial_survey_like(aerial_switch, aerial_years, aerial_obs, aerial_cv, aerial_cov, par_aerial_tau, par_log_aerial_sel, number_ysa, weight_fya)
  lp_aerial <- x$lp
  aerial_pred <- x$pred
  aerial_resid <- x$resid
  lp_aerial_tau <- x$lp_aerial_tau
  lp_troll <- get_troll_like(troll_switch, troll_years, troll_obs, troll_sd, par_troll_tau, number_ysa)
  # lp_tags <- get_tag_like(tag_switch, min_K + 1, n_K, n_T, n_I, n_J, 
  #                         first_yr, M_a, hrate_ysa,
  #                         par_hstar_i, tag_release_cta + 1, tag_recap_ctaa + 1,
  #                         minI = tag_rel_min_age + 1, 
  #                         maxI = tag_rel_max_age + 1, 
  #                         maxJ = tag_recap_max_age + 1,
  #                         shed1 = tag_shed_immediate, shed2 = tag_shed_continuous,
  #                         tag_rep_rates_ya,
  #                         tag_H_factor = par_tag_H_factor, tag_var_factor, tag_offset)
  # sum(lp_tags)
  # 176.553
  lp_tags <- 0
  # tag_pred <- array(0, dim = c(n_K, n_T, n_I, n_J))
  # tag_resid <- array(0, dim = c(n_K, n_T, n_I, n_J))
  lp_pop <- get_POP_like(pop_switch, pop_obs, phi_ya, spawning_biomass_y)
  lp_hsp <- get_HSP_like(hsp_switch, hsp_obs, par_hsp_q, hsp_false_negative, number_ysa, phi_ya, M_a, spawning_biomass_y, hrate_ysa)
  lp_gt <- get_GT_like(gt_switch, gt_obs, number_ysa)

  nll <- sum(lp_sel) + lp_rec + lp_hstar + lp_m10 + lp_h + lp_cpue_omega +
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
  REPORT(par_hstar_i)
  
  REPORT(lp_sel)
  REPORT(lp_rec)
  REPORT(lp_hstar)
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
  
  # REPORT(cpue_sigma)
  # REPORT(cpue_omega)
  REPORT(cpue_pred)
  REPORT(cpue_resid)
  REPORT(aerial_pred)
  REPORT(aerial_resid)
  # REPORT(troll_pred)
  # REPORT(troll_resid)
  # REPORT(tag_pred)
  # REPORT(tag_resid)
  # REPORT(lf_pred)
  # REPORT(af_pred)
  
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
  REPORT(tau_ac2)
  REPORT(recruitment_y)
  REPORT(par_rdev_y)
  # REPORT(rec_dev_y)
  REPORT(rdev_y)
  ADREPORT(par_sigma_r)

  return(nll)
}

#' Helper to make closure
#' 
#' @param f a \code{vector} of midpoints.
#' @param d natural mortality at the reference size.
#' @return a \code{vector}.
#' @export
#' 
cmb <- function(f, d) function(p) f(p, d)
