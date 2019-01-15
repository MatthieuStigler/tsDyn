library(tsDyn)
library(tidyverse)

data(barry)

## grid
grid_simple <- crossing(lag  =c(1, 2), 
                        include = c("const", "trend", "none", "both"))

## linear
models_VAR <-  grid_simple %>% 
  mutate(model = "VAR",
         object = map2(lag, include, ~ lineVar(barry, lag =.x, include = .y)))

models_VECM <-  grid_simple %>% 
  mutate(model = "VECM",
         object = map2(lag, include, ~ VECM(barry, lag =.x, include = .y)))

## setar
grid_setar <- crossing(lag  =c(1, 2), 
                       include = c("const", "trend", "none", "both"),
                       nthresh = 1:2)

models_TVAR <-  grid_setar %>% 
  mutate(model = "TVAR",
         object = pmap(list(lag, include, nthresh), 
                       ~suppressWarnings(TVAR(barry, lag =..1, include = ..2, nthresh=..3, 
                                               trace = FALSE))))

models_TVECM <-  grid_setar %>% 
  mutate(model = "TVECM",
         object = pmap(list(lag, include, nthresh), 
                       ~suppressWarnings(TVECM(barry[, 1:2], lag =..1, include = ..2, nthresh=..3, 
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
