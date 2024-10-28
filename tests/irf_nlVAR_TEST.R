library(tsDyn)
suppressMessages(library(dplyr))
library(purrr)
library(tidyr)
select <- dplyr::select
suppressWarnings(RNGversion("3.5.3"))

############################
### Load data
############################

path_mod_multi <- system.file("testdata/models_multivariate.rds", package = "tsDyn")

models_multivariate <- readRDS(path_mod_multi)

models_multivariate %>% 
  mutate(across(where(is.list), class)) %>% 
  as.data.frame() %>% 
  head(12)


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
  mutate(irf = map(object, ~irf_1(.)))

models_IRF_1$irf %>% 
  bind_rows() %>% 
  head() %>% 
  print(digits=3)

### irf_any
# irf.NULL <- function(x) NULL
# irf.ca.jo <- function(x) irf(vec2var(ca.jo))

models_VAR <- models_multivariate %>% 
  filter(model == "VAR")

## older method
models_IRF_any <- models_multivariate %>% 
  filter(model == "VAR") %>% 
  mutate(ortho = list(tibble(ortho =c(TRUE, FALSE)))) %>% 
  unnest(., ortho) %>% 
  mutate(irf = map2(object, ortho, ~irf_any(.x,  boot = TRUE, runs = 1, seed = 7, ortho = .y)),
         irf_vars = map2(object_vars, ortho, ~irf(.x, runs = 1, seed = 7, ortho = .y)),
         irf_vec2 = map2(object, ortho, ~irf(.x,  boot = FALSE, runs = 1, seed = 7, ortho = .y)))

models_IRF_any %>% 
  mutate(across(where(is.list), class)) %>% 
  as.data.frame()

## showquick summary
irf_extract_here <- function(x) {
  head(x$irf[[1]], 2) %>% 
    as.data.frame() %>% 
    mutate(type = "irf") %>% 
    rbind(head(x$Upper[[1]], 2) %>% 
            as.data.frame() %>% 
            mutate(type = "Upper_CI")) %>% 
    relocate(type)
}

## show head of irf any
map_dfr(models_IRF_any$irf, irf_extract_here) %>% 
  as.data.frame() %>% 
  head(10)%>% 
  mutate(across(where(is.numeric), ~round(., 6)))


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
         comp_irf_tsDOld_tsDNew = map2_lgl(irf, irf_vec2,  ~all.equal(.x$irf, .y$irf)),
         is_same_tsD_2ver = map_lgl(comp_irf_tsDOld_tsDNew, ~isTRUE(.))) %>% 
  dplyr::select(-starts_with("irf"), -starts_with("comp_irf"), comp_irf_tsDOld_tsDNew)

comp %>% 
  dplyr::select(-starts_with("object")) %>% 
  as.data.frame()

############################
### VECM
############################

models_VECM <- models_multivariate %>% 
  filter(model == "VECM") %>% 
  mutate(irf = map(object, ~irf_any(.,  boot = TRUE, runs = 1, seed = 7, ortho = FALSE)))

## show two first of first componment
models_VECM %>% 
  mutate(irf = map(irf, irf_extract_here)) %>% 
  dplyr::select(-object, -object_vars) %>% 
  unnest(irf) %>% 
  as.data.frame() %>% 
  mutate(across(where(is.numeric), ~round(., 6)))


## plot 1
plot(models_VECM$irf[[1]])
  


############################
### TVAR
############################


models_TVAR <- models_multivariate %>% 
  filter(model == "TVAR")

## test 1
tvar_1 <- models_TVAR$object[[1]]

irf(tvar_1, runs = 2, seed = 123)

## regime specific for TVAR
models_TVAR_irf <- models_TVAR  %>% 
  mutate(irf_L = map(object, ~irf_any(.,  boot = TRUE, runs = 1, seed = 7, ortho = FALSE, regime = "L")))

## show two first of first componment
models_TVAR_irf %>% 
  mutate(irf = map(irf_L, irf_extract_here)) %>%
  dplyr::select(-object, -object_vars, -irf_L ) %>% 
  unnest(irf) %>% 
  as.data.frame() %>% 
  mutate(across(where(is.numeric), ~round(., 6)))


## plot 1
plot(models_TVAR_irf$irf_L[[1]])

############################
### TVECM
############################

models_TVECM <- models_multivariate %>% 
  filter(model == "TVECM")

## test 1
tvecm_1 <- models_TVECM$object[[1]]
tsDyn:::irf_1(x=tvecm_1 , n.ahead = 10, cumulative = FALSE, regime = "L", ortho = TRUE)
tsDyn:::irf_1(x=tvecm_1 , n.ahead = 10, cumulative = FALSE, regime = "L", ortho = FALSE)
irf(x=tvecm_1, runs = 2, seed = 123)

## regime specific for TVECM
models_TVECM_irf <- models_TVECM   %>% 
  mutate(irf_L = map(object, ~suppressWarnings(irf_any(.,  boot = TRUE, runs = 1, seed = 7, ortho = FALSE, regime = "L"))))

## show two first of first componment
models_TVECM_irf %>% 
  mutate(irf = map(irf_L, irf_extract_here)) %>% 
  select(-object, -object_vars, -irf_L ) %>% 
  unnest(irf) %>% 
  as.data.frame() %>% 
  mutate(across(where(is.numeric), ~round(., 6)))


## plot 1
plot(models_TVECM_irf$irf_L[[1]])
