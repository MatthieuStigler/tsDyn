library(tsDyn)
suppressMessages(library(tidyverse))
suppressWarnings(RNGversion("3.5.3"))

############################
### Load data
############################
path_mod_multi <- system.file("inst/testdata/models_multivariate.rds", package = "tsDyn")
if(path_mod_multi=="") path_mod_multi <- system.file("testdata/models_multivariate.rds", package = "tsDyn")

models_multivariate <- readRDS(path_mod_multi)


mods <- models_multivariate$object
mods_nonLIn <- subset(models_multivariate, model %in% c("TVAR", "TVECM"))$object

############################
### tests
############################


## Standard functions
sapply(mods, class)
sapply(mods, print)
sapply(mods, summary)

sapply(mods, coef)
sapply(mods, tsDyn:::coefMat.nlVar)
sapply(mods, tsDyn:::coefVec.nlVar)
sapply(mods_nonLIn, coef, regime = "L")
sapply(mods_nonLIn, coef, regime = "H")

## utilities: get_lag, get_nVar, df.residual
sapply(mods, tsDyn:::get_series)
models_multivariate %>% 
  mutate(lag_out = map_int(object, tsDyn:::get_lag.nlVar),
         nVar_out = map_int(object, tsDyn:::get_nVar.nlVar),
         df_out = map(object, df.residual)) %>% 
  select(-starts_with("object")) %>% 
  as.data.frame()




uni_stats <- models_multivariate %>% 
  mutate_at("object", funs(deviance = map_dbl(., deviance),
                           AIC = map_dbl(., AIC),
                           BIC = map_dbl(., BIC),
                           logLik = map_dbl(., logLik))) %>% 
  select(-starts_with("object"))

as.data.frame(uni_stats)



sapply(mods, function(x) dim(residuals(x, initVal=FALSE)))
sapply(mods, function(x) dim(residuals(x, initVal=TRUE)))
sapply(mods, function(x) head(residuals(x), 3))
sapply(mods, function(x) head(residuals(x, initVal=TRUE), 3))
sapply(mods, function(x) tail(residuals(x), 3))

sapply(mods, function(x) head(fitted(x), 3))
sapply(mods, function(x) tail(fitted(x), 3))


##
suppressMessages(suppressWarnings(sapply(mods, tsDyn:::mod_refit_check)))

## NOn linear functions
sapply(mods_nonLIn, function(x) head(regime(x), 3))
sapply(mods_nonLIn, function(x) tail(regime(x), 3))

sapply(mods_nonLIn, function(x) head(regime(x, initVal=FALSE), 3))
sapply(mods_nonLIn, function(x) tail(regime(x, initVal=FALSE), 3))

sapply(mods_nonLIn, function(x) head(regime(x, time=FALSE), 3))
sapply(mods_nonLIn, function(x) head(regime(x, time=FALSE, initVal=FALSE), 3))


## toLatex
sapply(mods, toLatex)
options(show.signif.stars=FALSE)
sapply(mods, function(x) toLatex(summary(x)))

