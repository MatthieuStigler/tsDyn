library(tsDyn)
suppressMessages(library(tidyverse))

############################
### Load data
############################

if(getRversion()>="4.3.0"){
  file <- "models_multivariate.rds"
} else {
  file <- "models_multivariate_R_less_4_3.rds"
}

path_mod_multi <- system.file(file.path("inst/testdata/",file), package = "tsDyn")
if(path_mod_multi=="") path_mod_multi <- system.file(file.path("testdata/",file), package = "tsDyn")

models_multivariate <- readRDS(path_mod_multi)


mods <- models_multivariate$object
mods_nonLIn <- subset(models_multivariate, model %in% c("TVAR", "TVECM"))$object
mods_Linear <- subset(models_multivariate, model %in% c("VAR", "VECM"))$object
length(mods_nonLIn)+length(mods_Linear)==length(mods)

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
  mutate_at("object", list(deviance = ~map_dbl(., deviance),
                           AIC = ~map_dbl(., AIC),
                           BIC = ~map_dbl(., BIC),
                           logLik = ~map_dbl(., logLik))) %>% 
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

## Linear models
sapply(mods_Linear, \(x) round(confint(x),3))

## Non linear functions
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

