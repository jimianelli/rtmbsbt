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
  
  sel_min_age_f = c(2, 2, 2, 8, 6, 0, 2), # REMOVE f AS DIMENSION
  sel_max_age_f = c(17, 9, 17, 22, 25, 7, 17),
  sel_end_f = c(1, 0, 1, 1, 1, 0, 1),
  sel_change_sd_fy = t(as.matrix(sel_change_sd[,-1])), # CHANGE TO sel_change_year_fy
  
  pop_switch = 1, 
  hsp_switch = 1, hsp_false_negative = 0.7467647, 
  gt_switch = 1,
  cpue_switch = 1, cpue_a1 = 5, cpue_a2 = 17,
  aerial_switch = 4, aerial_tau = 0.3, 
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

min_len <- 86
bin_width <- 4
nbins <- 25
fsh <- data.frame(ifishery = 7, fishery = "CPUE", season = 2)

lf_data <- cpue_lfs <- read_csv("lf_assessment.csv") %>% filter(Fishery == 7)

obs_len_freq_il <- matrix(0, nrow = nrow(lf_data), ncol = 25)

for (irec in 1:nrow(lf_data)) {
  kbin <- 1
  mod_bin_wid <- min_len + bin_width * (kbin - 1)
  for (i in 1:110) {
    obs_bin_wid <- 32 + 2 * (i - 1);
    if (obs_bin_wid > mod_bin_wid && kbin < nbins) {
      kbin <- kbin + 1
      mod_bin_wid <- min_len + bin_width * (kbin - 1)
    }
    obs_len_freq_il[irec, kbin] <- obs_len_freq_il[irec, kbin] + as.numeric(lf_data[irec,-c(1:3)][i])
  }
  obs_len_freq_il[irec,] <- obs_len_freq_il[irec,] / sum(obs_len_freq_il[irec,])
}

ll <- seq(from = min_len, by = bin_width, length.out = nbins + 1) - 1

lf1 <- lf_data %>%
  pivot_longer(cols = !c(1:3), names_to = "Bin") %>%
  mutate(Bin = as.numeric(Bin)) %>%
  mutate(dBin = cut(Bin, breaks = ll, include.lowest = TRUE, right = FALSE)) %>%
  filter(!is.na(dBin)) %>%
  group_by(Fishery, Year, N, dBin) %>%
  summarise(value = sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = dBin) %>%
  left_join(fsh, by = join_by(Fishery == ifishery)) %>%
  relocate(Fishery, Year, season, N) %>%
  select(-fishery)

c(data$first_yr:data$last_yr)[data$lf_year[data$lf_fishery == 1]]
c(data$first_yr:data$last_yr)[data$cpue_years]
data$cpue_lfs <- obs_len_freq_il
data$cpue_n <- lf_data$N

source("../R/rtmb_functions.R")
xx <- get_selectivity(data$n_age, data$max_age, data$first_yr, data$first_yr_catch, 
                      data$sel_min_age_f, data$sel_max_age_f, data$sel_end_f, data$sel_change_year_fy,
                      data_par1$par_sels_init_i, data_par1$par_sels_change_i)
xx2 <- array(0, dim = c(7, 92, 31))
xx2[1:6,,] <- xx
xx2[7,,] <- xx[1,,]
xx <- xx2
data$sel_change_year_fy <- rbind(data$sel_change_year_fy, data$sel_change_year_fy[1,])
data$sel_change_year_fy[,data$first_yr_catch - data$first_yr + 1] <- 1 # must stay here for now
data$sel_change_year_fy[7,1:38] <- 0
data$sel_change_year_fy[7,]

par_sel <- list()
for (f in 1:7) {
  iy <- as.logical(data$sel_change_year_fy[f,])
  ia <- c(data$sel_min_age_f[f]:data$sel_max_age_f[f]) + 1
  par_sel[[f]] <- log(xx[f, iy, ia])
}
par_sel[[4]] <- t(as.matrix(par_sel[[4]])) # to force as matrix
plot(xx[1,22,], type = "l")
for (i in 55:75) lines(xx[1,i,], col = i)
plot(exp(par_sel[[1]][1,]), type = "l")
for (i in 2:10) lines(exp(par_sel[[1]][i,]), col = i)

parameters <- list(
  par_log_B0 = 16.19836, par_log_psi = log(1.5),
  par_log_m0 = log(0.4), par_log_m4 = log(0.1670507),
  par_log_m10 = log(0.065), par_log_m30 = log(0.45741),
  par_log_h = log(0.55),
  # par_log_h = log(0.8),
  par_log_sigma_r = log(0.6), 
  par_log_cpue_q = -0.02033773, par_cpue_creep = 0.005,
  par_log_cpue_sigma = log(0.2), par_log_cpue_omega = log(1),
  par_log_aerial_tau = log(0.3), par_log_aerial_sel = c(0, 0),
  par_log_troll_tau = log(0.3689704), par_log_hsp_q = 0, 
  par_log_tag_H_factor = log(1),
  par_sel_rho_y = c(0.7, 0.7, 0.5, 0.7, 0.5, 0.5, 0.7),
  par_sel_rho_a = c(0.9, 0.9, 0.5, 0.9, 0.9, 0.5, 0.9),
  par_log_sel_sigma = log(c(0.31, 0.25, 0.75, 0.25, 0.38, 1.13, 0.31)),
  par_log_sel_1 = par_sel[[1]], par_log_sel_2 = par_sel[[2]],
  par_log_sel_3 = par_sel[[3]], par_log_sel_4 = par_sel[[4]],
  par_log_sel_5 = par_sel[[5]], par_log_sel_6 = par_sel[[6]],
  par_log_sel_7 = par_sel[[7]], 
  par_rdev_y = data_par1$Reps
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
map[["par_sel_rho_y"]] <- factor(rep(NA, 7))
map[["par_sel_rho_a"]] <- factor(rep(NA, 7))
map[["par_log_sel_sigma"]] <- factor(rep(NA, 7))

source("../R/model.R")
source("../R/rtmb_functions.R")
obj <- RTMB::MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)
# obj <- RTMB::MakeADFun(func = cmb(sbt_model, data), 
#                        parameters = parameters, map = map, random = c("par_log_sel_6"))
unique(names(obj$par))
obj$report()$lp_lf
obj$fn()

Params <- parameters
bnd <- get_bounds(obj = obj) # NEEDS WORK AND THE parameters AS INPUT

control <- list(eval.max = 10000, iter.max = 10000)

# opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
#               lower = bnd$lower, upper = bnd$upper, control = control)
opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he,
              lower = bnd$lower, upper = bnd$upper, control = control)
opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he,
              lower = bnd$lower, upper = bnd$upper, control = control)

ce <- check_estimability(obj = obj)
ce[[4]] %>% filter(Param_check != "OK")

# Check OK to use 2DAR1 for single year ----

f <- 4 # ADD TO TEST
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
p1 <- plot_selectivity(data = data, object = obj, years = 1969:2022, fisheries = "CPUE")
p2 <- plot_selectivity(data = data, object = obj, years = 1969:2022, fisheries = "LL1")
p1 + p2
plot_cpue_lf(data = data, object = obj)

plot_selectivity(data = data, object = obj)
plot_selectivity(data = data, object = obj, years = 1954:1991, fisheries = "LL4")

plot_lf(data = data, object = obj, fishery = "LL1")
plot_lf(data = data, object = obj, fishery = "LL4")
plot_af(data = data, object = obj, fishery = "Indonesian")
plot_af(data = data, object = obj, fishery = "Australian")
plot_cpue(data = data, object = obj)
plot_biomass_spawning(data = data, object = obj)

plot(obj$report()$number_ysa[39,2,], type = "l")
plot(obj$report()$sel_fya[7,39,], type = "l")
plot(obj$report()$number_ysa[39,2,] * obj$report()$sel_fya[7,39,], type = "l")
plot(obj$report()$cpue_lf_pred[1,], type = "l")
plot(obj$report()$spawning_biomass_y, type = "l")

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
                 get_recruitment = get_recruitment, 
                 get_harvest_rate = get_harvest_rate,
                 get_rho = get_rho, get_selectivity2 = get_selectivity2,
                 get_sel_like = get_sel_like, get_tag_like = get_tag_like,
                 get_recruitment_prior = get_recruitment_prior,
                 get_length_like = get_length_like, get_age_like = get_age_like, 
                 get_cpue_like = get_cpue_like, get_troll_like = get_troll_like,
                 get_cpue_length_like = get_cpue_length_like, 
                 get_aerial_survey_like = get_aerial_survey_like, 
                 get_POP_like = get_POP_like, get_HSP_like = get_HSP_like, 
                 get_GT_like = get_GT_like))

save(data, parameters, obj, opt, mcmc, file = "mcmc_0divergences.rda")
plot_sampler_params(fit = mcmc, plot = TRUE)
decamod::pairs_rtmb(fit = mcmc, order = "slow", pars = 1:5)
decamod::pairs_rtmb(fit = mcmc, order = "mismatch", pars = 1:5)
decamod::pairs_rtmb(fit = mcmc, order = "fast", pars = 1:5)
library(ellipse)
library(GGally)
source("../../../decamod/R/plots.R")
pairs_rtmb(fit = mcmc, order = "divergent", pars = 1:5)
plot_uncertainties(fit = mcmc, log = TRUE, plot = TRUE)


##--Jim's tests--------------
##
##


case=1
sel_phi <- matrix(0, nrow = 6, ncol = 2)
sel_phi[,1] <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.5) # year
sel_phi[,2] <- c(0.9, 0.9, 0.9, 0.9, 0.9, 0.5) # age
sel_scale <- matrix(0, nrow = 6, ncol = 2)
sel_scale[,1] <- c(0.8, 0.8, 0.8, 0.8, 1.2, 2) # year
sel_scale[,2] <- c(0.8, 0.8, 0.8, 0.8, 1.2, 2) # year
df <- list()
df[[1]] <- scale_phi(case = case, sel_phi, sel_scale) 

case=2
sel_phi <- matrix(0, nrow = 6, ncol = 2)
sel_scale <- matrix(0, nrow = 6, ncol = 2)
sel_phi[,1] <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.5) # year
sel_phi[,2] <- c(0.9, 0.9, 0.9, 0.9, 0.9, 0.5) # age
sel_scale[,1] <- c(1, 0.8, 1.0, 0.8, 1.2, 1.5) # year
sel_scale[,2] <- c(1, 0.8, 1.0, 0.8, 1.2, 1.5) # age   sel_scale[,2] <- c(0.8, 0.8, 0.8, 0.8, 0.5, 2) # age
df[[case]] <- scale_phi(case = case, sel_phi, sel_scale) 

case=3
sel_phi <- matrix(0, nrow = 6, ncol = 2)
sel_scale <- matrix(0, nrow = 6, ncol = 2)
sel_phi[,1] <- c(0.5, 0.7, 0.5, 0.7, 0.7, 0.5) # year
sel_phi[,2] <- c(0.5, 0.9, 0.5, 0.9, 0.9, 0.5) # age
sel_scale[,1] <- c(1, 0.8, 1.0, 0.8, 1.2, 1.5) # year
sel_scale[,2] <- c(1, 0.8, 1.0, 0.8, 1.2, 1.5) # age   sel_scale[,2] <- c(0.8, 0.8, 0.8, 0.8, 0.5, 2) # age
df[[case]] <- scale_phi(case = case, sel_phi, sel_scale)      

case=4 # Case 2 ll1, 3 for aus and LL3, try to get Indonesian better
sel_phi <- matrix(0, nrow = 6, ncol = 2)
sel_phi[,1] <- c(0.7, 0.7, 0.5, 0.7, 0.5, 0.5) # year
sel_phi[,2] <- c(0.9, 0.9, 0.5, 0.9, 0.9, 0.5) # age
sel_scale <- matrix(0, nrow = 6, ncol = 2)
sel_scale[,1] <- c(1, 0.8, 1.0, 0.8, 1.2, 1.5) # year
sel_scale[,2] <- c(1, 0.8, 1.0, 0.8, 1.2, 1.5) # age   sel_scale[,2] <- c(0.8, 0.8, 0.8, 0.8, 0.5, 2) # age
df[[case]] <- scale_phi(case = case, sel_phi, sel_scale)      

case=5 # Case 2 ll1, 3 for aus and LL3, try to get Indonesian better
sel_phi <- matrix(0, nrow = 6, ncol = 2)
sel_scale <- matrix(0, nrow = 6, ncol = 2)
sel_phi[,1] <- c(0.7, 0.7, 0.5, 0.7, 0.5, 0.5) # year
sel_phi[,2] <- c(0.9, 0.9, 0.5, 0.9, 0.9, 0.5) # age
sel_scale[,1] <- c(1, 0.8, 1.0, 0.8, 1.0, 1.5) # year
sel_scale[,2] <- c(1, 0.8, 1.0, 0.8, 1.0, 1.5) # age   sel_scale[,2] <- c(0.8, 0.8, 0.8, 0.8, 0.5, 2) # age
df[[case]] <- scale_phi(case = case, sel_phi, sel_scale)      

bind_rows(df) |> ggplot(aes(x = fishery, y = RTMB_nll, fill = as.factor(case))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(ADMB_nll, 0)), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Comparison of RTMB and ADMB NLL by Fishery and Case",
       x = "Fishery", y = "Negative Log-Likelihood (NLL)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, size=.5, hjust = 1))
names(bind_rows(df))
bind_rows(df) |> select(ggplot(aes(x = fishery, y = RTMB_nll, fill = as.factor(case))) +

sel_scale
scale_phi <- function(case,sel_phi=sel_phi,sel_scale=sel_scale) {
  parameters$par_log_sel_phi   <- log(sel_phi)
  parameters$par_log_sel_scale <- log(sel_scale)
  obj <- RTMB::MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)
  unique(names(obj$par))
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he,
                lower = bnd$lower, upper = bnd$upper, control = control)
  opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he,
                lower = bnd$lower, upper = bnd$upper, control = control)
  obj$par <- opt$par
  
  af_df <- data.frame(fishery = data$af_fishery, value = obj$report()$lp_af) %>%
    group_by(fishery) %>%
    summarise(RTMB_nll = sum(value), .groups = "drop") %>%
    bind_cols(
      case = case,
      ADMB_nll = lr$lnlike[5:6],
      RTMB_sel = obj$report()$lp_sel[5:6],
      phi_yr = sel_phi[5:6, 1],
      phi_age = sel_phi[5:6, 2],
      scale_yr = sel_scale[5:6, 1],
      scale_age = sel_scale[5:6, 2]
    )
  # Second block (lf), renamed to match
  lf_df <- data.frame(fishery = data$lf_fishery, value = obj$report()$lp_lf) %>%
    group_by(fishery) %>%
    summarise(RTMB_nll = sum(value), .groups = "drop") %>%
    bind_cols(
      case = case,
      ADMB_nll = lr$lnlike[1:4],
      RTMB_sel = obj$report()$lp_sel[1:4],
      phi_yr = sel_phi[1:4, 1],
      phi_age = sel_phi[1:4, 2],
      scale_yr = sel_scale[1:4, 1],
      scale_age = sel_scale[1:4, 2]
    )
  # Combine
  phi_scale_df <- bind_rows(af_df, lf_df)
  phi_scale_df
}

