
library(RTMB)
library(TMBhelper)
library(reshape2)
library(tidyverse)
library(sbt)
library(kableExtra)
library(testthat)

theme_set(theme_bw())

attach(data_csv1)
lr <- data_labrep1

data <- list(
  last_yr = 2022, age_increase_M = 25,
  length_m50 = 150, length_m95 = 180, 
  catch_UR_on = 0, catch_surf_case = 1, catch_LL1_case = 1, 
  scenarios_surf = scenarios_surface, scenarios_LL1 = scenarios_LL1,
  sel_min_age_f = c(2, 2, 2, 8, 6, 0), 
  sel_max_age_f = c(17, 9, 17, 22, 25, 7),
  sel_end_f = c(1, 0, 1, 1, 1, 0),
  sel_change_sd_fy = t(as.matrix(sel_change_sd[,-1])), 
  sel_smooth_sd_f = lr$sel.smooth.sd,
  pop_switch = 1, 
  hsp_switch = 1, hsp_false_negative = 0.7467647, 
  gt_switch = 1,
  cpue_switch = 1, cpue_a1 = 5, cpue_a2 = 17,
  aerial_switch = 4, aerial_tau = data_labrep1$tau.aerial, 
  troll_switch = 1, 
  lf_minbin = c(1, 1, 1, 11),
  tag_switch = 1, tag_var_factor = 1.82
)

data <- get_data(data_in = data)

data$cpue_years <- data$cpue_years + 1
data$aerial_years <- data$aerial_years + 1
data$troll_years <- data$troll_years + 1
data$af_year <- data$af_year + 1
data$lf_year <- data$lf_year + 1
data$pop_obs[,1] <- data$pop_obs[,1] + 1
data$sel_change_year_fy[,1] <- 1

n_sely_f <- rowSums(data$sel_change_year_fy)
n_sely_f
n_sela_f <- sel_max_age_f - sel_min_age_f + 1
n_sela_f
par_sel <- list()
for (f in 1:6) par_sel[[f]] <- matrix(0, nrow = n_sely_f[f], ncol = n_sela_f[f])

parameters <- list(
  par_log_B0 = data_par1$ln_B0,
  par_log_psi = log(data_par1$psi),
  par_log_m0 = log(data_par1$m0), 
  par_log_m4 = log(data_par1$m4),
  par_log_m10 = log(data_par1$m10), 
  par_log_m30 = log(data_par1$m30),
  par_log_h = log(data_par1$steep), 
  par_log_sigma_r = log(lr$sigma.r), 
  par_rdev_y = data_par1$Reps,
  # par_sels_init_i = data_par1$par_sels_init_i, 
  # par_sels_change_i = data_par1$par_sels_change_i,
  par_log_sel_1 = par_sel[[1]],
  par_log_sel_2 = par_sel[[2]],
  par_log_sel_3 = par_sel[[3]],
  par_log_sel_4 = par_sel[[4]],
  par_log_sel_5 = par_sel[[5]],
  par_log_sel_6 = par_sel[[6]],
  par_log_cpue_q = data_par1$lnq,
  par_log_cpue_sigma = log(data_par1$sigma_cpue),
  par_log_cpue_omega = log(data_par1$cpue_omega),
  par_log_aerial_tau = log(data_par1$tau_aerial),
  par_log_aerial_sel = data_par1$ln_sel_aerial,
  par_log_troll_tau = log(data_par1$tau_troll),
  par_log_hsp_q = data_par1$lnqhsp, 
  par_logit_hstar_i = qlogis(exp(data_par1$par_log_hstar_i)),
  par_log_tag_H_factor = log(data_par1$tag_H_factor)
)

map <- list()
map[["par_log_psi"]] <- factor(NA)
map[["par_log_m0"]] <- factor(NA)
map[["par_log_m10"]] <- factor(NA)
map[["par_log_h"]] <- factor(NA)
map[["par_log_sigma_r"]] <- factor(NA)
map[["par_log_cpue_sigma"]] <- factor(NA)
map[["par_log_cpue_omega"]] <- factor(NA)
map[["par_log_aerial_tau"]] <- factor(NA)
map[["par_log_aerial_sel"]] <- factor(rep(NA, 2))
map[["par_log_hsp_q"]] <- factor(NA)
map[["par_log_tag_H_factor"]] <- factor(NA)

Params <- parameters
bnd <- get_bounds(obj = obj)

source("../R/model.R")
source("../R/rtmb_functions.R")
obj <- RTMB::MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)
# note that when catch and rec devs are all set to zero the SSB is not flat.

control <- list(eval.max = 10000, iter.max = 10000)

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he,
              lower = bnd$lower, upper = bnd$upper, control = control)
opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he,
              lower = bnd$lower, upper = bnd$upper, control = control)

opt$convergence
ce <- check_estimability(obj = obj)
ce[[4]] %>% filter(Param_check != "OK")
# obj$par[1:10]
# obj$report()$spawning_biomass_y
# opt$par[1:10]
# obj$par <- opt$par

library(adnuts)
mcmc <- sample_sparse_tmb(
  obj = obj, metric = "auto", iter = 1000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.95), init = "last.par.best",
  # lower = Lwr, upper = Upr, # these bounds dont seem to work
  globals = list(posfun = posfun, get_M = get_M, get_phi = get_phi, 
                 get_initial_numbers = get_initial_numbers, 
                 get_recruitment = get_recruitment, get_harvest_rate = get_harvest_rate,
                 get_rho = get_rho, get_selectivity = get_selectivity,
                 get_sel_like = get_sel_like, get_recruitment_prior = get_recruitment_prior,
                 get_length_like = get_length_like, get_age_like = get_age_like, 
                 get_cpue_like = get_cpue_like, get_tag_like = get_tag_like,
                 get_aerial_survey_like = get_aerial_survey_like, get_troll_like = get_troll_like,
                 get_POP_like = get_POP_like, get_HSP_like = get_HSP_like, get_GT_like = get_GT_like))

# save(mcmc, file = "mcmc_no1.rda")
plot_sampler_params(fit = mcmc, plot = TRUE)
decamod::pairs_rtmb(fit = mcmc, order = "slow", pars = 1:5)
decamod::pairs_rtmb(fit = mcmc, order = "mismatch", pars = 1:5)
decamod::pairs_rtmb(fit = mcmc, order = "fast", pars = 1:5)
plot_uncertainties(fit = mcmc, log = TRUE, plot = TRUE)
