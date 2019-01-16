library(tsDyn)
library(tidyverse)

############################
### Load data
############################
path_mod_multi <- system.file("inst/testdata/models_multivariate.rds", package = "tsDyn")
if(path_mod_multi=="") path_mod_multi <- system.file("testdata/models_multivariate.rds", package = "tsDyn")

models_multivariate <- readRDS(path_mod_multi)

models_multivariate

############################
### VAR
############################

irf_any <- tsDyn:::irf_any
irf_1 <- tsDyn:::irf_1
irf_1.nlVar <- tsDyn:::irf_1.nlVar

## manual comparisons
mod_random_1 <- filter(models_multivariate, lag ==2)$object[[2]]
mod_random_1_vars <- filter(models_multivariate, lag ==2)$object_vars[[2]]

irf_any(mod_random_1, boot = FALSE)$irf[[1]]
irf(mod_random_1, boot = FALSE)$irf[[1]]
irf(mod_random_1_vars, boot = FALSE)$irf[[1]]

irf_any(mod_random_1, boot = FALSE, ortho = FALSE)$irf[[1]]
irf(mod_random_1, boot = FALSE, ortho = FALSE)$irf[[1]]
irf(mod_random_1_vars, boot = FALSE, ortho = FALSE)$irf[[1]]

### irf _1
models_IRF_1 <- models_multivariate %>% 
  filter(model == "VAR") %>% 
  mutate(irf = map(object, ~irf_1(.,  boot = TRUE, runs = 2, seed = 7)))

models_IRF_1$irf %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  head()

### irf_any
# irf.NULL <- function(x) NULL
# irf.ca.jo <- function(x) irf(vec2var(ca.jo))

models_VAR <- models_multivariate %>% 
  filter(model == "VAR")

## older method
models_IRF_any <- models_multivariate %>% 
  filter(model == "VAR") %>% 
  mutate(ortho = list(tibble(ortho =c(TRUE, FALSE)))) %>% 
  unnest(ortho, .preserve = c("object", "object_vars")) %>% 
  mutate(irf = map2(object, ortho, ~irf_any(.x,  boot = TRUE, runs = 2, seed = 7, ortho = .y)),
         irf_vars = map2(object_vars, ortho, ~irf(.x, runs = 2, seed = 7, ortho = .y)),
         irf_vec2 = map2(object, ortho, ~irf(.x,  boot = FALSE, runs = 2, seed = 7, ortho = .y)))

models_IRF_any

## show head of irf any
map_df(models_IRF_any$irf, ~ head(.$irf[[1]], 2) %>%  as_tibble) %>% 
  as.data.frame() %>% 
  as_tibble()


## compare with vars
all.equal(models_IRF_any$irf[[1]]$irf, 
          models_IRF_any$irf_vars[[1]]$irf)
models_IRF_any$irf[[1]]$irf[[1]]
models_IRF_any$irf_vars[[1]]$irf[[1]]
models_IRF_any$irf_vec2[[1]]$irf[[1]]

comp <- models_IRF_any %>% 
  mutate(comp_irf_tsD_vars = map2(irf, irf_vars,  ~all.equal(.x$irf, .y$irf)),
         is_same = map_lgl(comp_irf_tsD_vars, ~isTRUE(.)),
         comp_irf_tsDOld_vars = map2(irf_vec2, irf_vars,  ~all.equal(.x$irf, .y$irf)),
         is_same_tssDvec2 = map_lgl(comp_irf_tsDOld_vars, ~isTRUE(.)),
         comp_irf_tsDOld_tsDNew = map2(irf, irf_vec2,  ~all.equal(.x$irf, .y$irf)),
         is_same_tsD_2ver = map_lgl(comp_irf_tsDOld_tsDNew, ~isTRUE(.))) %>% 
  select(-starts_with("irf"), -starts_with("comp_irf"), comp_irf_tsDOld_tsDNew)

comp %>% 
  select(-starts_with("object"))

## regime specific for TVAR
models_TVAR <- models_multivariate %>% 
  filter(model == "TVAR") 

models_TVAR %>% 
  mutate(irf_L = map(object, ~irf_any(.,  boot = TRUE, runs = 2, seed = 7, ortho = FALSE, regime = "L")))
