library(tsDyn)
library(vars)
library(tidyverse)

data(barry)

## grid
grid_simple <- crossing(lag  =c(1, 2), 
                        include = c("const", "trend", "none", "both"))

data_inp <-  barry[, 1:2]

## VAR, VECM
var_upd <- function(data, p, type){
  res <- VAR(data, p = p, type = type)
  res$call$type <- type
  res$call$p <- p
  res
}

models_VAR <-  grid_simple %>% 
  mutate(model = "VAR",
         object = map2(lag, include, ~ lineVar(data_inp, lag =.x, include = .y)),
         object_vars = map2(lag, include, ~ var_upd(data_inp, .x, .y)))

models_VAR$object_vars[[1]]$call


models_VECM_tsD <-  grid_simple %>% 
  mutate(model = "VECM",
         object = map2(lag, include, ~ VECM(data_inp, lag =.x, include = .y)))

models_VECM_vars <- grid_simple %>% 
  filter(include == "const") %>% 
  mutate(model = "VECM",
         object_vars = map2(lag, include, ~ ca.jo(data_inp, K =.x+1, spec = "transitory")))

models_VECM <- models_VECM_tsD %>% 
  left_join(models_VECM_vars, by = c("lag", "include", "model"))

## TVAR, TVEC<
grid_tvar <- crossing(lag  =c(1, 2), 
                       include = c("const", "trend", "none", "both"),
                       nthresh = 1:2)

models_TVAR <-  grid_tvar %>% 
  mutate(model = "TVAR",
         object = pmap(list(lag, include, nthresh), 
                       ~suppressWarnings(TVAR(data_inp, lag =..1, include = ..2, nthresh=..3, 
                                               trace = FALSE))))

models_TVECM <-  grid_tvar %>% 
  mutate(model = "TVECM",
         object = pmap(list(lag, include, nthresh), 
                       ~suppressWarnings(TVECM(data_inp, lag =..1, include = ..2, nthresh=..3, 
                                              trace = FALSE, plot = FALSE))))



## combine all

models_multivariate <- bind_rows(models_VAR, 
                                 models_VECM,
                                 models_TVAR,
                                 models_TVECM)

## saving it in sysdata, see http://r-pkgs.had.co.nz/data.html#data-sysdata
## no!! would be loaded
## so save in: inst/testdata https://stackoverflow.com/questions/32328802/where-to-put-data-for-automated-tests-with-testthat
path <- system.file("inst/testdata", package = "tsDyn")
if(path== "") path <- system.file("testdata", package = "tsDyn")
path <-  "~/Dropbox/Documents/tsDyn/tsDyn/inst/testdata"
saveRDS(models_multivariate, file= paste(path, "models_multivariate.rds", sep = "/"))

## this gives path: system.file("inst/testdata/models_multivariate.rds", package = "tsDyn")
path_mod_multi <- system.file("inst/testdata/models_multivariate.rds", package = "tsDyn")
if(path_mod_multi=="") path_mod_multi <- system.file("testdata/models_multivariate.rds", package = "tsDyn")
