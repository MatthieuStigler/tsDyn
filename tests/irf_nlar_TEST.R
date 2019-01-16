library(tsDyn)
suppressMessages(library(tidyverse))


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
df_regs <-  tibble(model = c("linear", "setar", "setar"),
                       regime = c("all", "L", "H"))
models_irf <- models_univariate %>% 
  filter(!model %in% c("aar", "lstar" )) %>% 
  left_join(df_regs, by = "model") %>% 
  mutate(irf = map2(object, regime,  ~suppressWarnings(irf(.x,  boot = TRUE, runs = 20, seed = 7, regime = .y))))

## IRF
df_irf <- map_df(models_irf$irf, ~ head(.$irf[[1]], 2) %>%  as_tibble) %>% 
  as.data.frame()

## Lower
df_all <- models_irf %>% 
  mutate(irf_irf = map(irf, ~ head(.$irf[[1]], 5)),
         irf_low = map(irf, ~ head(.$Lower[[1]], 5)),
         irf_upp = map(irf, ~ head(.$Upper[[1]], 5))) %>% 
  select(-irf) %>% 
  gather(irf_stat, value, irf_irf, irf_low, irf_upp) %>% 
  mutate(value = map(value, ~tibble(x=.) %>% 
                          mutate(n.ahead = 0:4))) %>% 
  unnest(value) %>% 
  spread(irf_stat, x)

df_all %>% 
  filter(n.ahead %in% c( 1)) %>% 
  as.data.frame()


df_all %>% 
  mutate(is_in = irf_irf >= irf_low & irf_irf <= irf_upp) %>% 
  count(model, regime, is_in)

