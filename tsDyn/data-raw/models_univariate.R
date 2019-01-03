library(tsDyn)
library(tidyverse)

## grid
grid_simple <- crossing(lag  =c(1, 2), 
                        include = c("const", "trend", "none", "both"))

## linear
models_linear <-  grid_simple %>% 
  mutate(model = "linear",
         object = map2(lag, include, ~ linear(lh, m =.x, include = .y)))

## setar
grid_setar <- crossing(lag  =c(1, 2), 
                       include = c("const", "trend", "none", "both"),
                       nthresh = 1:2, 
                       thDelay = 0:1)

models_setar <-  grid_setar %>% 
  filter(thDelay<lag) %>% 
  mutate(model = "setar",
         object = pmap(list(lag, include, nthresh, thDelay), 
                       ~suppressWarnings(setar(lh, m =..1, include = ..2, nthresh=..3, 
                                               thDelay = ..4,
                                               trace = FALSE))))


## lstar
grid_lstar <- crossing(lag  =c(1, 2), 
                       include = c("const", "trend", "none", "both"),
                       thDelay = 0:1)

models_lstar <-  grid_lstar %>% 
  filter(thDelay<lag) %>% 
  mutate(model = "lstar",
         object = pmap(list(lag, include, thDelay), 
                       ~suppressWarnings(lstar(lh, m =..1, include = ..2, thDelay=..3, trace = FALSE))))

## aar
grid_aar <- crossing(lag  =c(1, 2, 3))

models_aar <-  grid_aar %>% 
  mutate(model = "aar",
         object = map(lag, ~suppressWarnings(aar(lh, m =.))))

## combine all

models_univariate <- bind_rows(models_linear, 
                               models_setar,
                               models_lstar,
                               models_aar)

## saving it in sysdata, see http://r-pkgs.had.co.nz/data.html#data-sysdata
## no!! would be loaded
## so save in: inst/testdata https://stackoverflow.com/questions/32328802/where-to-put-data-for-automated-tests-with-testthat
path <- paste(system.file("inst/testdata", package = "tsDyn"), "models_univariate.rds", sep = "/")
saveRDS(models_univariate, file=path)

## this gives path: system.file("inst/testdata/models_univariate.rds", package = "tsDyn")
path_mod_uni <- system.file("inst/testdata/models_univariate.rds", package = "tsDyn")
if(path_mod_uni=="") path_mod_uni <- system.file("testdata/models_univariate.rds", package = "tsDyn")
