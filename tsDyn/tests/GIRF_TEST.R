
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

models_ar_setar


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


