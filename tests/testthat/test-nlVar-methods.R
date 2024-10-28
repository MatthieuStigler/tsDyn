library(tsDyn)
suppressMessages(library(dplyr))
library(purrr)

############################
### Load data
############################

if(getRversion()>="4.3.0"){
  file <- "models_multivariate.rds"
} else {
  file <- "models_multivariate_R_less_4_3.rds"
}

# path_mod_multi <- testthat::test_path("testdata", file)
path_mod_multi <- system.file("testdata/models_multivariate.rds", package = "tsDyn")
models_multivariate <- readRDS(path_mod_multi)


mods <- models_multivariate$object
mods_nonLIn <- subset(models_multivariate, model %in% c("TVAR", "TVECM"))$object
mods_Linear <- subset(models_multivariate, model %in% c("VAR", "VECM"))$object
length(mods_nonLIn)+length(mods_Linear)==length(mods)


################################
#'## tests
################################

## Standard functions
test_that("simple methods", {
  expect_snapshot(sapply(mods, class))
  expect_snapshot(sapply(mods, print))
  expect_snapshot(sapply(mods, summary))
  
  expect_snapshot(sapply(mods, coef))
  expect_snapshot(sapply(mods, tsDyn:::coefMat.nlVar))
  expect_snapshot(sapply(mods, tsDyn:::coefVec.nlVar))
  expect_snapshot(sapply(mods_nonLIn, coef, regime = "L"))
  expect_snapshot(sapply(mods_nonLIn, coef, regime = "H"))
})


## 
test_that("utilities: get_lag, get_nVar, df.residual", {
  expect_snapshot(sapply(mods, tsDyn:::get_series))
  expect_snapshot(models_multivariate %>% 
    mutate(lag_out = map_int(object, tsDyn:::get_lag.nlVar),
           nVar_out = map_int(object, tsDyn:::get_nVar.nlVar),
           df_out = map(object, df.residual)) %>% 
    select(-starts_with("object")) %>% 
    as.data.frame())
})


test_that("AIC et al", {
  uni_stats <- models_multivariate %>% 
    mutate_at("object", list(deviance = ~map_dbl(., deviance),
                                             AIC = ~map_dbl(., AIC),
                                             BIC = ~map_dbl(., BIC),
                                             logLik = ~map_dbl(., logLik))) %>% 
                      select(-starts_with("object"))
  
  expect_snapshot(as.data.frame(uni_stats))
})

test_that("More methods", {
  
  expect_snapshot(sapply(mods, function(x) dim(residuals(x, initVal=FALSE))))
  expect_snapshot(sapply(mods, function(x) dim(residuals(x, initVal=TRUE))))
  expect_snapshot(sapply(mods, function(x) head(residuals(x), 3)))
  expect_snapshot(sapply(mods, function(x) head(residuals(x, initVal=TRUE), 3)))
  expect_snapshot(sapply(mods, function(x) tail(residuals(x), 3)))
  
  expect_snapshot(sapply(mods, function(x) head(fitted(x), 3)))
  expect_snapshot(sapply(mods, function(x) tail(fitted(x), 3)))
})

##
test_that("mod_refit_check", {
  expect_snapshot(suppressMessages(suppressWarnings(sapply(mods, tsDyn:::mod_refit_check))))
})


## Linear models
test_that("Linear models", {
  expect_snapshot(sapply(mods_Linear, \(x) round(confint(x),3)))
  expect_snapshot(map_dfr(mods_Linear, \(x) tidy(x) |> 
                            mutate(across(where(is.numeric), ~round(., 3)))))
})

## Non linear functions
test_that("Non linear functions", {
  expect_snapshot(sapply(mods_nonLIn, function(x) head(regime(x), 3)))
  expect_snapshot(sapply(mods_nonLIn, function(x) tail(regime(x), 3)))
  
  expect_snapshot(sapply(mods_nonLIn, function(x) head(regime(x, initVal=FALSE), 3)))
  expect_snapshot(sapply(mods_nonLIn, function(x) tail(regime(x, initVal=FALSE), 3)))
  
  expect_snapshot(sapply(mods_nonLIn, function(x) head(regime(x, time=FALSE), 3)))
  expect_snapshot(sapply(mods_nonLIn, function(x) head(regime(x, time=FALSE, initVal=FALSE), 3)))
})

## toLatex
test_that("toLatex", {
  expect_snapshot(sapply(mods, toLatex))
  expect_snapshot(options(show.signif.stars=FALSE))
  expect_snapshot(sapply(mods, function(x) toLatex(summary(x))))
})
