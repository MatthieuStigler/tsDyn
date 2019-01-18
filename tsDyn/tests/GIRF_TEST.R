
library(tsDyn)
suppressMessages(library(tidyverse))

############################
### Load data
############################
path_mod_uni <- system.file("inst/testdata/models_univariate.rds", package = "tsDyn")
if(path_mod_uni=="") path_mod_uni <- system.file("testdata/models_univariate.rds", package = "tsDyn")

path_mod_multi <- system.file("inst/testdata/models_multivariate.rds", package = "tsDyn")
if(path_mod_multi=="") path_mod_multi <- system.file("testdata/models_multivariate.rds", package = "tsDyn")

models_ar_setar <- readRDS(path_mod_uni) %>% 
  filter(model %in% c("ar", "setar"))

models_multivariate <- readRDS(path_mod_multi)

############################
### Univariate
############################


mod_1_uni <- models_ar_setar$object[[1]]

tsDyn:::irf_1_shock(mod_1_uni, 
                    shock = 1,
                    hist = 0,
                    seed = 123)

tsDyn:::irf_1_shock_ave(object = mod_1_uni, 
                        shock = 1,
                        hist = 0, 
                        seed = 123)


GIRF(object = mod_1_uni, 
     hist_li = list(rep(1.6, 1)),
     shock_li = list(0.01), 
     R = 20, 
     seed = 123) %>% head


## Simple, given shocks
models_ar_setar %>% 
  head(2) %>% 
  mutate(girf = map2(object, lag, ~GIRF(object=.x, n.ahead = 3,
                                        hist_li = list(rep(1.6, .y)),
                                        shock_li = list(0.01), R = 20, seed = 123) %>% as_tibble)) %>% 
  unnest(girf) %>% 
  as.data.frame()

## Simple, random
models_ar_setar %>% 
  head(2) %>% 
  mutate(girf = map(object, ~GIRF(object=., n.ahead = 3, n.hist = 3, n.shock = 3,
                                  R = 20, seed = 123) %>% as_tibble)) %>% 
  unnest(girf) %>% 
  as.data.frame()


############################
### Multivariate
############################

mod_1_multi <- models_multivariate$object[[1]]

tsDyn:::irf_1_shock(mod_1_multi, 
                    shock = matrix(c(1, 0), nrow = 1),
                    hist = matrix(c(0, 0), nrow = 1))

tsDyn:::irf_1_shock_ave(mod_1_multi, 
                        shock = matrix(c(1, 0), nrow = 1),
                        hist = matrix(c(0, 0), nrow = 1),
                        seed = 123)

GIRF(object=mod_1_multi, 
     shock_li = list(matrix(c(1, 0), nrow = 1),
                     matrix(c(0, 1), nrow = 1)),
     hist_li = list(matrix(c(0, 0), nrow = 1),
                    matrix(c(0, 1), nrow = 1)),
     seed = 123) %>% 
  head
