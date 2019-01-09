library(tsDyn)
library(tidyverse)

############################
### Load data
############################
path_mod_uni <- system.file("inst/testdata/models_univariate.rds", package = "tsDyn")
if(path_mod_uni=="") path_mod_uni <- system.file("testdata/models_univariate.rds", package = "tsDyn")

models_univariate <- readRDS(path_mod_uni)

############################
### Compute
############################

## ar_mean
models_ar_mean <- models_univariate %>% 
  filter(model!="aar") %>% 
  mutate(ar_mean = map(object, ~suppressWarnings(ar_mean(.))))

############################
### Show result
############################

## should be only for both and trend
check_1 <- models_ar_mean %>% 
  filter(map_lgl(ar_mean, is.null)) %>% 
  count(include)

stopifnot(all(check_1$include == c("both", "trend")))

## should be only for const and none
check_2 <- models_ar_mean %>% 
  filter(!map_lgl(ar_mean, is.null)) %>% 
  count(include)

stopifnot(all(check_2$include == c("const", "none")))

## show results
models_ar_mean %>% 
  filter(!map_lgl(ar_mean, is.null)) %>% 
  mutate(ar_mean = map(ar_mean, ~as_tibble(t(.)))) %>% 
  unnest(ar_mean) %>% 
  select(-object) %>% 
  print(n = Inf)

