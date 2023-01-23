library(tsDyn)
suppressWarnings(RNGversion("3.5.3"))


############################
### Load data
############################
path_mod_uni <- system.file("inst/testdata/models_univariate.rds", package = "tsDyn")
if(path_mod_uni=="") path_mod_uni <- system.file("testdata/models_univariate.rds", package = "tsDyn")

models_univariate <- readRDS(path_mod_uni)



mod <-  models_univariate$object
mod_no_aar <-  subset(models_univariate,  model != "aar")$object
mod_notrend_noaar <-  subset(models_univariate,  !include %in% c("trend", "both") & model != "aar")$object
mod_notrend <-  subset(models_univariate,  !include %in% c("trend", "both") )$object
mod_const_only <-  subset(models_univariate,  include =="const" )$object

### Extract methods
sapply(mod, AIC)
sapply(mod, BIC)
sapply(mod, mse)
sapply(mod, MAPE)

sapply(mod, coef)
sapply(mod, coef, hyperCoef = FALSE)
sapply(mod, function(x) head(residuals(x)))
sapply(mod, function(x) head(residuals(x, initVal = FALSE)))


lapply(mod_const_only, predict, n.ahead=10)

lapply(mod_const_only, predict, n.ahead=3, type="MC", seed=1234)
lapply(mod_const_only, predict, n.ahead=3, type="bootstrap", seed=1234)
lapply(mod_const_only, predict, n.ahead=3, type="block-bootstrap", seed=1234)

## charac root
lapply(mod_notrend_noaar, charac_root)
lapply(mod_notrend_noaar, ar_mean)

### Utility functions
sapply(mod, getTh)


## Output of mod_no_aar[-44] is platform/machien specific...
## Output of mod_no_aar[-23] is platform/machien specific: doesn't work on M1mac
suppressMessages(suppressWarnings(sapply(mod_no_aar[-c(23, 44)], tsDyn:::mod_refit_check)))


### Pred Roll, acc_stat:
x <- log10(lynx)
mod <- list()
mod[["linear"]] <- linear(x, m=2)
mod[["setar"]] <- setar(x, m=2, thDelay=1, trace = FALSE)
mod[["star"]] <- star(x, m=2, thDelay=1, trace = FALSE)
mod[["lstar"]] <- lstar(x, m=2, thDelay=1, trace = FALSE)
mod[["aar"]] <- aar(x, m=2)


x_small <- x[1:100]
mod_small <- list()
mod_small[["linear"]] <- linear(x_small, m=2)
mod_small[["setar"]] <- setar(x_small, m=2, thDelay=1, th=getTh(mod[["setar"]]), trace = FALSE)
mod_small[["lstar"]] <- lstar(x_small, m=2, thDelay=1, th=getTh(mod[["lstar"]]), gamma=coef(mod[["lstar"]])["gamma"], trace = FALSE)
mod_small[["aar"]] <- aar(x_small, m=2)

pred_rolls_1 <- lapply(mod_small, predict_rolling, n.ahead=1, newdata=x[101:114])
sapply(pred_rolls_1, function(x) x$pred[[1]])
sapply(pred_rolls_1, accuracy_stat)[-1,] ## removing first line as gave 'factor' under R<4


pred_rolls_12 <- lapply(mod_small, predict_rolling, n.ahead=1:2, newdata=x[101:114])
sapply(pred_rolls_12, function(x) x$pred[[1]])
lapply(pred_rolls_12, accuracy_stat)

