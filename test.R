library(rtmbFunctions)

# Example usage of rho_trans
rho_trans(0.5)

# Example usage of double_logistic
x <- seq(0, 10, length.out = 100)
y <- double_logistic(x, p1 = 1, p2 = 2, p3 = 3)
plot(x, y, type = 'l')

# Example usage of get_sel_skeleton
age_a <- 1:10
par_log_sel_skel <- matrix(rnorm(18), nrow = 6, ncol = 3)
sel_skeleton_fa <- get_sel_skeleton(age_a, par_log_sel_skel)

# Example usage of get_selectivity
sel_switch_f <- c(2, 2, 2, 2, 2, 2)
n_year <- 10
n_age <- 10
max_age <- 10
first_yr <- 1931
first_yr_catch <- 1952
first_yr_catch_f <- c(1952, 1952, 1952, 1952, 1952, 1952)
sel_min_age_f <- c(1, 1, 1, 1, 1, 1)
sel_max_age_f <- c(10, 10, 10, 10, 10, 10)
sel_end_f <- c(1, 1, 1, 1, 1, 1)
sel_change_year_fy <- matrix(0, nrow = 6, ncol = 10)
sel_change_year_fy[, 5] <- 1
par_log_sel1_ay <- array(rnorm(100), dim = c(10, 10))
par_log_sel2_ay <- array(rnorm(100), dim = c(10, 10))
par_log_sel3_ay <- array(rnorm(100), dim = c(10, 10))
par_log_sel4_ay <- array(rnorm(100), dim = c(10, 10))
par_log_sel5_ay <- array(rnorm(100), dim = c(10, 10))
par_log_sel6_ay <- array(rnorm(100), dim = c(10, 10))
sel_fya <- get_selectivity(sel_switch_f, n_year, n_age, max_age, age_a, first_yr, first_yr_catch, 
                           first_yr_catch_f, sel_min_age_f, sel_max_age_f, sel_end_f, 
                           sel_change_year_fy, sel_skeleton_fa, par_log_sel1_ay, par_log_sel2_ay, 
                           par_log_sel3_ay, par_log_sel4_ay, par_log_sel5_ay, par_log_sel6_ay)

# Example usage of construct_Q
n_years <- 10
n_ages <- 10
ay_Index <- matrix(1:(n_years * n_ages), nrow = n_years * n_ages, ncol = 2)
rho_y <- 0.5
rho_a <- 0.5
rho_c <- 0.5
log_sigma2 <- 0.1
Var_Param <- 0
Q_sparse <- construct_Q(n_years, n_ages, ay_Index, rho_y, rho_a, rho_c, log_sigma2, Var_Param)
warnings()
