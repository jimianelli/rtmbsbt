library(RTMB)
library(TMBhelper)
library(reshape2)
library(tidyverse)
library(sbt)
library(kableExtra)
library(testthat)

source("../R/get-data.R")

theme_set(theme_bw())

attach(data_csv1)
lr <- data_labrep1

data <- list(
  last_yr = 2022, age_increase_M = 25,
  length_m50 = 150, length_m95 = 180, 
  catch_UR_on = 0, catch_surf_case = 1, catch_LL1_case = 1, 
  scenarios_surf = scenarios_surface, scenarios_LL1 = scenarios_LL1,
  removal_switch_f = c(0, 0, 0, 1, 0, 0), # 0=harvest rate, 1=direct removals
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
# data$pop_obs[,1] <- data$pop_obs[,1] + 1
# data$sel_max_age_f[4] <- 18

source("../R/rtmb_functions.R")
xx <- get_selectivity(data$n_age, data$max_age, data$first_yr, data$first_yr_catch, 
                      data$sel_min_age_f, data$sel_max_age_f, data$sel_end_f, data$sel_change_year_fy,
                      data_par1$par_sels_init_i, data_par1$par_sels_change_i)

data$sel_change_year_fy[,data$first_yr_catch - data$first_yr + 1] <- 1 # must stay here for now

par_sel <- list()
for (f in 1:6) {
  iy <- as.logical(data$sel_change_year_fy[f,])
  ia <- c(data$sel_min_age_f[f]:data$sel_max_age_f[f]) + 1
  par_sel[[f]] <- log(xx[f,iy,ia])
}
par_sel[[4]] <- t(as.matrix(par_sel[[4]]))
plot(xx[1,22,], type = "l")
for (i in 55:75) lines(xx[1,i,], col = i)

plot(exp(par_sel[[1]][1,]), type = "l")
for (i in 2:10) lines(exp(par_sel[[1]][i,]), col = i)

# sel_phi <- matrix(0, nrow = 6, ncol = 2)
# sel_scale <- matrix(0, nrow = 6, ncol = 2)
# sel_phi[,1] <- c(0.7, 0.7, 0.5, 0.7, 0.5, 0.5) # year
# sel_phi[,2] <- c(0.9, 0.9, 0.5, 0.9, 0.9, 0.5) # age
# sel_scale[,1] <- c(1, 0.8, 1.0, 0.8, 1.0, 1.5) # year
# sel_scale[,2] <- c(1, 0.8, 1.0, 0.8, 1.0, 1.5) # age
rho_y <- c(0.7, 0.7, 0.5, 0.7, 0.5, 0.5) # year
rho_a <- c(0.9, 0.9, 0.5, 0.9, 0.9, 0.5) # age
sigma <- c(1, 0.8, 1.0, 0.8, 1.0, 1.5) * sqrt(1 - rho_y^2) * sqrt(1 - rho_a^2)
(scale <- sqrt(sigma^2) / sqrt(1 - rho_y^2) / sqrt(1 - rho_a^2))

parameters <- list(
  par_log_B0 = data_par1$ln_B0,
  par_log_psi = log(data_par1$psi),
  par_log_m0 = log(data_par1$m0), 
  par_log_m4 = log(data_par1$m4),
  par_log_m10 = log(data_par1$m10), 
  par_log_m30 = log(data_par1$m30),
  par_log_h = log(data_par1$steep),
  # par_log_h = log(0.8),
  par_log_sigma_r = log(lr$sigma.r), 
  par_rdev_y = data_par1$Reps,
  par_log_sel_1 = par_sel[[1]],
  par_log_sel_2 = par_sel[[2]],
  par_log_sel_3 = par_sel[[3]],
  par_log_sel_4 = par_sel[[4]],
  par_log_sel_5 = par_sel[[5]],
  par_log_sel_6 = par_sel[[6]],
  par_log_sel_7 = par_sel[[7]],
  par_sel_rho_y = rho_y,
  par_sel_rho_a = rho_a,
  par_log_sel_sigma = log(sigma),
  # par_log_sel_phi = log(sel_phi),
  # par_log_sel_scale = log(sel_scale),
  par_log_cpue_q = data_par1$lnq,
  par_cpue_creep = 0.005,
  par_log_cpue_sigma = log(data_par1$sigma_cpue),
  par_log_cpue_omega = log(data_par1$cpue_omega),
  par_log_aerial_tau = log(data_par1$tau_aerial),
  par_log_aerial_sel = data_par1$ln_sel_aerial,
  par_log_troll_tau = log(data_par1$tau_troll),
  par_log_hsp_q = data_par1$lnqhsp, 
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
map[["par_cpue_creep"]] <- factor(NA)
map[["par_log_aerial_tau"]] <- factor(NA)
map[["par_log_aerial_sel"]] <- factor(rep(NA, 2))
map[["par_log_hsp_q"]] <- factor(NA)
map[["par_log_tag_H_factor"]] <- factor(NA)
map[["par_log_sel_4"]] <- factor(matrix(NA, nrow = nrow(parameters$par_log_sel_4), ncol = ncol(parameters$par_log_sel_4)))
map[["par_sel_rho_y"]] <- factor(rep(NA, 6))
map[["par_sel_rho_a"]] <- factor(rep(NA, 6))
map[["par_log_sel_sigma"]] <- factor(rep(NA, 6))
# map_phi <- matrix(NA, nrow = 6, ncol = 2)
# # map_phi[6,] <- c(1, 2)
# map[["par_log_sel_phi"]] <- factor(map_phi)
# map_scale <- matrix(NA, nrow = 6, ncol = 2)
# # map_scale[6,] <- c(1, 1)
# map[["par_log_sel_scale"]] <- factor(map_scale)

source("../R/model.R")
source("../R/rtmb_functions.R")
obj <- RTMB::MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)
# obj <- RTMB::MakeADFun(func = cmb(sbt_model, data), 
#                        parameters = parameters, map = map, random = c("par_log_sel_6"))
unique(names(obj$par))
obj$report()$lp_pop
obj$fn()

Params <- parameters
bnd <- get_bounds(obj = obj)

control <- list(eval.max = 10000, iter.max = 10000)

# opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
#               lower = bnd$lower, upper = bnd$upper, control = control)
opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he,
              lower = bnd$lower, upper = bnd$upper, control = control)
opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he,
              lower = bnd$lower, upper = bnd$upper, control = control)

ce <- check_estimability(obj = obj)
ce[[4]] %>% filter(Param_check != "OK")
# obj$par[1:10]
# obj$report()$spawning_biomass_y
# opt$par[1:10]
# obj$par <- opt$par

# Check OK to use 2DAR1 for single year ----

f <- 4
sigma2 <- exp(parameters$par_log_sel_sigma[f])^2
scale <- sqrt(sigma2) / sqrt(1 - rho_y^2) / sqrt(1 - rho_a^2) # Define 2d scale
f1 <- function(x) dautoreg(x, phi = 0.7, log = TRUE) # year
f2 <- function(x) dautoreg(x, phi = 0.9, log = TRUE) # age
dseparable(f1, f2)(par_sel[[4]], scale = scale)
dautoreg(par_sel[[4]][1,], phi = 0.9, log = TRUE, scale = scale)

exp(obj$par[names(obj$par) %in% c("par_log_sel_phi", "par_log_sel_scale")])

source("../../sbt/R/plot-selectivity.R")
library(scales)
library(ggridges)
plot_selectivity(data = data, object = obj)
plot_selectivity(data = data, object = obj, years = 1954:1991, fisheries = "LL4")
plot_lf(data = data, object = obj, fishery = "LL1")
plot_lf(data = data, object = obj, fishery = "LL4")
plot_af(data = data, object = obj, fishery = "Indonesian")
plot_af(data = data, object = obj, fishery = "Australian")

data.frame(fishery = data$af_fishery, value = obj$report()$lp_af) %>% 
  group_by(fishery) %>% 
  summarise(RTMB_nll = sum(value)) %>%
  bind_cols(ADMB_nll = lr$lnlike[5:6], RTMB_sel = obj$report()$lp_sel[5:6])

data.frame(fishery = data$lf_fishery, value = obj$report()$lp_lf) %>% 
  group_by(fishery) %>% 
  summarise(RTMB = sum(value)) %>%
  bind_cols(ADMB = lr$lnlike[1:4])

library(adnuts)
mcmc <- sample_sparse_tmb(
  obj = obj, metric = "auto", iter = 1000, warmup = 750, chains = 4, cores = 4,
  # obj = obj, metric = "auto", iter = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.99), init = "last.par.best",
  # lower = Lwr, upper = Upr, # these bounds dont seem to work
  globals = list(posfun = posfun, get_M = get_M, get_phi = get_phi, 
                 get_initial_numbers = get_initial_numbers, 
                 get_recruitment = get_recruitment, get_harvest_rate = get_harvest_rate,
                 get_rho = get_rho, 
                 get_selectivity = get_selectivity, get_selectivity2 = get_selectivity2,
                 get_sel_like = get_sel_like, get_recruitment_prior = get_recruitment_prior,
                 get_length_like = get_length_like, get_age_like = get_age_like, 
                 get_cpue_like = get_cpue_like, get_tag_like = get_tag_like,
                 get_aerial_survey_like = get_aerial_survey_like, get_troll_like = get_troll_like,
                 get_POP_like = get_POP_like, get_HSP_like = get_HSP_like, get_GT_like = get_GT_like))

save(data, parameters, obj, opt, mcmc, file = "mcmc__3divergences.rda")
plot_sampler_params(fit = mcmc, plot = TRUE)
decamod::pairs_rtmb(fit = mcmc, order = "slow", pars = 1:5)
decamod::pairs_rtmb(fit = mcmc, order = "mismatch", pars = 1:5)
decamod::pairs_rtmb(fit = mcmc, order = "fast", pars = 1:5)
library(ellipse)
library(GGally)
source("../../../decamod/R/plots.R")
pairs_rtmb(fit = mcmc, order = "divergent", pars = 1:5)
plot_uncertainties(fit = mcmc, log = TRUE, plot = TRUE)
