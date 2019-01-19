
library(tsDyn)
suppressMessages(library(tidyverse))


plot_GIRF_line_low <-  tsDyn:::plot_GIRF_line_low
irf_1_shock <- tsDyn:::irf_1_shock
irf_1_shock_ave <- tsDyn:::irf_1_shock_ave



add.regime = FALSE
shock_both = TRUE
n.hist =  5
n.ahead =  10
n.shock = 20

############################
### Load data
############################
path_mod_uni <- system.file("inst/testdata/models_univariate.rds", package = "tsDyn")
if(path_mod_uni=="") path_mod_uni <- system.file("testdata/models_univariate.rds", package = "tsDyn")

path_mod_multi <- system.file("inst/testdata/models_multivariate.rds", package = "tsDyn")
if(path_mod_multi=="") path_mod_multi <- system.file("testdata/models_multivariate.rds", package = "tsDyn")

models_ar_setar <- readRDS(path_mod_uni) %>% 
  filter(model %in% c("linear", "setar"))

models_multivariate <- readRDS(path_mod_multi)

############################
### Univariate
############################


mod_1_uni <- models_ar_setar %>% 
  filter(model == "setar") %>% {.$object[[1]]}

mod_1_ar <- models_ar_setar %>% 
  filter(model == "linear") %>% {.$object[[1]]}


mod_1_uni_1_shock <- irf_1_shock(object=mod_1_uni, 
                                 shock = 1,
                                 hist = 0,
                                 seed = 123)

mod_ar_1_shock <- irf_1_shock(mod_1_ar, shock = 1, hist = 0,seed = 123)

mod_1_uni_1_shock
plot_GIRF_line_low(x=mod_1_uni_1_shock)

plot_GIRF_line_low(x=mod_ar_1_shock)
plot(irf(mod_1_ar, boot = FALSE))

mod_1_uni_1_shock_ave <- irf_1_shock_ave(object = mod_1_uni, 
                                         shock = 1,
                                         hist = 0, 
                                         seed = 123)

plot_GIRF_line_low(x=mod_1_uni_1_shock_ave)

mod_1_uni_1_girf <- GIRF(object = mod_1_uni, 
                         hist_li = list(1),
                         shock_li = list(0.1, 0.11, -0.1, -0.11), 
                         R = 2, 
                         seed = 123)

head(mod_1_uni_1_girf)
plot_GIRF_line_low(mod_1_uni_1_girf, n_simu  = 1:4)

mod_1_uni_1_girf_big <- GIRF(object = mod_1_uni, 
                             n.hist = 5,
                             R = 2, 
                             seed = 123)

plot(x=mod_1_uni_1_girf_big, plot_type = "density")
plot(x=mod_1_uni_1_girf_big, plot_type = "line", n_simu = 1:50, add_legend =FALSE)

## Simple, given shocks
models_ar_setar %>% 
  mutate(girf = map2(object, lag, ~GIRF(object=.x, n.ahead = 3,
                                        hist_li = list(rep(1.6, .y)),
                                        shock_li = list(0.01), R = 2, seed = 123) %>% as_tibble)) %>% 
  unnest(girf) %>% 
  as.data.frame()

## Simple, random
models_ar_setar %>% 
  mutate(girf = map(object, ~GIRF(object=., n.ahead = 3, n.hist = 3, n.shock = 3,
                                  R = 2, seed = 123) %>% as_tibble)) %>% 
  unnest(girf) %>% 
  as.data.frame()


############################
### Multivariate
############################

mod_TVAR <- models_multivariate %>% 
  filter(model == "TVAR" & lag ==2) %>% {.$object[[1]]}

mod_VAR <- models_multivariate %>% 
  filter(model == "VAR"& lag ==2) %>% {.$object[[1]]}

mod_VECM <- models_multivariate %>% 
  filter(model == "VECM" & lag ==2) %>% {.$object[[1]]}

mod_TVAR_1_shock <- irf_1_shock(object = mod_TVAR, 
                                   shock = matrix(c(1, 0), nrow = 1),
                                   hist = matrix(c(0, 0, 0, 0), nrow = 2),
                                   seed = 123)

mod_VAR_1_shock <- irf_1_shock(object = mod_VAR, shock = matrix(c(1, 0), nrow = 1),
                               hist = matrix(c(0, 0, 0, 0), nrow = 2), seed = 123)

mod_VECM_1_shock <- irf_1_shock(object = mod_VECM, shock = matrix(c(1, 0), nrow = 1),
                                hist = matrix(rep(0, 6), nrow = 3), seed = 123)


plot_GIRF_line_low(x=mod_VAR_1_shock)
plot_GIRF_line_low(x=mod_TVAR_1_shock)

plot_GIRF_line_low(x=mod_TVAR_1_shock, var = "cpiUSA")
plot_GIRF_line_low(x=mod_VAR_1_shock, var = "cpiUSA")

mod_TVAR_1_shock_ave <- irf_1_shock_ave(mod_TVAR, 
                                        shock = matrix(c(1, 0), nrow = 1),
                                        hist = matrix(c(0, 0, 0, 0), nrow = 2),
                                        seed = 123)


plot_GIRF_line_low(x=mod_TVAR_1_shock_ave)

TVAR_GIRF <- GIRF(object=mod_TVAR, 
                  shock_li = list(matrix(c(1, 0), nrow = 1),
                                  matrix(c(0, 1), nrow = 1)),
                  hist_li = list(matrix(c(0, 0, 0, 0), nrow = 2),
                                 matrix(c(0, 1, 0, 0), nrow = 2)),
                  R = 2,
                  seed = 123) 

gi_out <- GIRF(object=mod_TVAR, seed = 123, n.hist = 40, R = 2) 
plot(density(residuals(mod_TVAR)[, 1]))
head(gi_out)
plot(x=gi_out, var = "dolcan")
plot(x=gi_out, var = "cpiUSA")

## Simple, random
models_multivariate %>% 
  # head(2) %>% 
  mutate(girf = map(object, ~GIRF(object=., n.ahead = 3, n.hist = 3, n.shock = 3,
                                  R = 2, seed = 123) %>% head(2))) %>% 
  unnest(girf) %>% 
  as.data.frame()

