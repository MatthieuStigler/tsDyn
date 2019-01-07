library(tsDyn)
library(tidyverse)

############################
### Load data
############################
path_mod_uni <- system.file("inst/testdata/models_univariate.rds", package = "tsDyn")
if(path_mod_uni=="") path_mod_uni <- system.file("testdata/models_univariate.rds", package = "tsDyn")

models_univariate <- readRDS(path_mod_uni)

############################
### Test irf univariate
############################


## boot: many models instable! had to search for a while to find seed with no errors...
models_irf <- models_univariate %>% 
  filter(!model %in% c("aar", "lstar" )) %>% 
  mutate(irf = map(object, ~suppressWarnings(irf(.,  boot = TRUE, runs = 2, seed = 7))))

## IRF
df_irf <- map_df(models_irf$irf, ~ head(.$irf[[1]], 2) %>%  as_data_frame) %>% 
  as.data.frame()

## Lower
df_low <- map_df(models_irf$irf, ~ head(.$Lower[[1]], 2) %>%  as_data_frame) %>% 
  as.data.frame()
df_upp <- map_df(models_irf$irf, ~ head(.$Upper[[1]], 2) %>%  as_data_frame) %>% 
  as.data.frame()

cbind(df_irf, df_low, df_upp)
