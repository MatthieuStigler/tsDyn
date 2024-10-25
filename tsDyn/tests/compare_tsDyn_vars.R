library(broom)
library(tsDyn)
library(vars)
suppressMessages(library(dplyr))
library(purrr)

data(Canada)

options(useFancyQuotes=FALSE) # useful for all.equal comparison

#########################
##### VECM #####
#########################

### unrestricted cons
vecm_l1_co_tsD <-VECM(Canada, lag=1, include="const", estim="ML")
vecm_l3_co_tsD <-VECM(Canada, lag=3, include="const", estim="ML")

vecm_l1_co_var <- ca.jo(Canada, K=2, ecdet="none", spec="transitory")
vecm_l3_co_var <- ca.jo(Canada, K=4, ecdet="none", spec="transitory")

### unrestricted cons and trend
vecm_l1_bo_tsD <-VECM(Canada, lag=1, include="both", estim="ML")
vecm_l3_bo_tsD <-VECM(Canada, lag=3, include="both", estim="ML")

trend_it <- function(start = 3) matrix(c(rep(0, start-1), seq_len(nrow(Canada)-start+1)), ncol =1, dimnames = list(NULL, "trend"))
vecm_l1_bo_var <- ca.jo(Canada, K=2, ecdet="none", spec="transitory", dumvar = trend_it(3))
vecm_l3_bo_var <- ca.jo(Canada, K=4, ecdet="none", spec="transitory", dumvar = trend_it(5))


### restricted cons
vecm_l1_LRco_tsD <-VECM(Canada, lag=1, LRinclude="const", estim="ML", include="none")
vecm_l1_LRco_var <- ca.jo(Canada, K=2, ecdet="const", spec="transitory")

vecm_l3_LRco_tsD <-VECM(Canada, lag=3, LRinclude="const", estim="ML", include="none")
vecm_l3_LRco_var <- ca.jo(Canada, K=4, ecdet="const", spec="transitory")

### restricted trend
vecm_l1_LRtr_tsD <-VECM(Canada, lag=1, LRinclude="trend", estim="ML")
vecm_l1_LRtr_var <- ca.jo(Canada, K=2, ecdet="trend", spec="transitory")

vecm_l3_LRtr_tsD <-VECM(Canada, lag=3, LRinclude="trend", estim="ML")
vecm_l3_LRtr_var <- ca.jo(Canada, K=4, ecdet="trend", spec="transitory")

all_models <- list(
  l1_co = list(vecm_l1_co_var, vecm_l1_co_tsD), 
  l3_co = list(vecm_l3_co_var, vecm_l3_co_tsD), 
  l1_bo = list(vecm_l1_bo_var, vecm_l1_bo_tsD),
  l3_bo = list(vecm_l3_bo_var, vecm_l3_bo_tsD),
  l1_LRco = list(vecm_l1_LRco_var, vecm_l1_LRco_tsD),
  l3_LRco = list(vecm_l3_LRco_var, vecm_l3_LRco_tsD),
  l1_LRtr = list(vecm_l1_LRtr_var, vecm_l1_LRtr_tsD),
  l3_LRtr = list(vecm_l3_LRtr_var, vecm_l3_LRtr_tsD))

deftol <- .Machine$double.eps ^ 0.5
lowtol <- 1.e-07
lowlowtol <- 1.e-06
comp_teststat <- function(x, tol=deftol) all.equal(x[[1]]@teststat, rev(rank.test(x[[2]])$res_df[,"eigen"]), check.attributes=FALSE, tolerance=tol)
comp_betas <- function(x, tol=deftol) all.equal(cajorls(x[[1]])$beta, x[[2]]$model.specific$beta, check.attributes=FALSE, tolerance=tol)
comp_coefs <- function(x, tol=deftol) all.equal(coefficients(cajorls(x[[1]])$rlm), t(coefficients(x[[2]])), check.attributes=FALSE, tolerance=tol)
comp_LL <- function(x) all.equal(as.numeric(logLik(vec2var(x[[1]]))), logLik(x[[2]]), check.attributes=FALSE)
comp_IRF <- function(x, ortho = TRUE) all.equal(irf(vec2var(x[[1]]), boot=FALSE, ortho = ortho)$irf,
                                                irf(x[[2]], boot=FALSE, ortho = ortho)$irf, check.attributes=FALSE)
comp_FEVD <- function(x) all.equal(fevd(vec2var(x[[1]])), fevd(x[[2]]), check.attributes=FALSE)
comp_resid <- function(x, tol=deftol) all.equal(residuals(vec2var(x[[1]])), residuals(x[[2]]), check.attributes=FALSE, tolerance=tol)
comp_fitted <- function(x) all.equal(fitted(vec2var(x[[1]])), fitted(x[[2]], level="original"), check.attributes=FALSE)
comp_predictOld <- function(x) all.equal(predict(vec2var(x[[1]]))$fcst, tsDyn:::predictOld.VECM(x[[2]])$fcst, check.attributes=FALSE)
comp_predict <- function(x) all.equal(sapply(predict(vec2var(x[[1]]), n.ahead=5)$fcst,function(x) x[,"fcst"]), predict(x[[2]]), check.attributes=FALSE)
comp_VECM_coefA <- function(x,...) all.equal(coefA(x[[1]]), coefA(x[[2]]), check.attributes=FALSE,...)
comp_VECM_coefB <- function(x,...) all.equal(coefB(x[[1]]), coefB(x[[2]]), check.attributes=FALSE,...)
comp_VECM_coefPI <- function(x, ...) all.equal(coefPI(x[[1]]), coefPI(x[[2]]), check.attributes=FALSE, ...)
comp_VECM_logLik <- function(x, r=1) all.equal(logLik(x[[1]], r=r), logLik(x[[2]], r=r), check.attributes=FALSE)

### Small function to print nicely output of all.equal, rounding the number:

roundAll.Equal <- function(x, round=8){
  isFALSE <- x!="TRUE"
  xFalse <- x[isFALSE]
  # extract the number (i.e remoe all the rest)
  xf<- gsub("(Component ([0-9]+)?([[:punct:]][[:alnum:]]+[[:punct:]])?: )?Mean relative difference: ", 
            "", xFalse)
  xf2<- round(as.numeric(xf),round)
  x[isFALSE] <- paste("Mean relative difference: ", xf2, sep="")
  x
}

models_noLR <-  all_models[-grep("LR", names(all_models))]

### Compare VECM methods:
sapply(all_models, comp_teststat )
sapply(all_models, comp_betas, tol=lowtol)
roundAll.Equal(sapply(all_models, comp_coefs, tol = 1e-07), round=7) # 5 and 6
sapply(all_models, comp_LL)
sapply(models_noLR, comp_IRF)
sapply(models_noLR, comp_IRF, ortho = FALSE)
sapply(all_models, comp_FEVD)
sapply(all_models, comp_resid, tol=lowlowtol) # 5
sapply(all_models, comp_fitted)
sapply(all_models[-c(3, 4)], comp_predict)
lapply(sapply(all_models[-c(3, 4)], comp_predictOld),roundAll.Equal, round=6) # 5 and 6
sapply(all_models, comp_VECM_coefA, tolerance=1e-07)
sapply(all_models, comp_VECM_coefB, tolerance=1e-07)
sapply(all_models, comp_VECM_coefPI, tolerance=1e-07)
mapply(function(r) sapply(all_models, comp_VECM_logLik, r=r), r=0:4)

## restricted betas:
all.equal(coefA(vecm_l1_co_tsD),coefA(vecm_l1_co_var), check.attributes=FALSE)
all.equal(coefB(vecm_l1_co_tsD),coefB(vecm_l1_co_var), check.attributes=FALSE)

R <- matrix(c(1,0.1, -0.24, 3.6), ncol=1)

vecm_Rrest_tsD <-VECM(Canada, lag=1, include="const", estim="ML", beta=R)
vecm_Rrest_var <-blrtest(vecm_l1_co_var, H=R, r=1)

all.equal(coefA(vecm_Rrest_tsD),coefA(vecm_Rrest_var), check.attributes=FALSE)
all.equal(coefB(vecm_Rrest_tsD),coefB(vecm_Rrest_var), check.attributes=FALSE)
all.equal(logLik(vecm_Rrest_tsD, r=1),logLik(vecm_Rrest_var, r=1), check.attributes=FALSE)

all.equal(deviance(vecm_Rrest_tsD),sum(deviance(cajorls(vecm_Rrest_var)[[1]])), check.attributes=FALSE)
all.equal(-2*(logLik(vecm_Rrest_tsD)-logLik(vecm_l1_co_tsD)), vecm_Rrest_var@teststat)

#########################
##### VAR #####
#########################

var_l1_co_tsD <-lineVar(Canada, lag=1, include="const")
var_l1_tr_tsD <-lineVar(Canada, lag=1, include="trend")
var_l1_bo_tsD <-lineVar(Canada, lag=1, include="both")
var_l1_no_tsD <-lineVar(Canada, lag=1, include="none")

var_l3_co_tsD <-lineVar(Canada, lag=3, include="const")
var_l3_tr_tsD <-lineVar(Canada, lag=3, include="trend")
var_l3_bo_tsD <-lineVar(Canada, lag=3, include="both")
var_l3_no_tsD <-lineVar(Canada, lag=3, include="none")

var_l1_co_var <- VAR(Canada, p=1, type="const")
var_l1_tr_var <- VAR(Canada, p=1, type="trend")
var_l1_bo_var <- VAR(Canada, p=1, type="both")
var_l1_no_var <- VAR(Canada, p=1, type="none")

var_l3_co_var <- VAR(Canada, p=3, type="const")
var_l3_tr_var <- VAR(Canada, p=3, type="trend")
var_l3_bo_var <- VAR(Canada, p=3, type="both")
var_l3_no_var <- VAR(Canada, p=3, type="none")



all_var_models <- list(
  l1_co = list(var_l1_co_tsD, var_l1_co_var), 
  l1_tr= list(var_l1_tr_tsD, var_l1_tr_var), 
  l1_bo= list(var_l1_bo_tsD, var_l1_bo_var), 
  l1_no= list(var_l1_no_tsD, var_l1_no_var), 
  l3_co= list(var_l3_co_tsD, var_l3_co_var), 
  l3_tr= list(var_l3_tr_tsD, var_l3_tr_var), 
  l3_bo= list(var_l3_bo_tsD, var_l3_bo_var), 
  l3_no= list(var_l3_no_tsD, var_l3_no_var))

all_var_models_noNoBo <- all_var_models[-grep("bo|no", names(all_var_models))]


## test functions
coef_to_vars <- function(x){
  c <- coef(x)
  if(any(grepl("Intercept|Trend", colnames(c)))){
    wh.deter <- grep("Intercept|Trend", colnames(c))
    res <- cbind(c[,-wh.deter], c[,wh.deter,drop=FALSE])
  } else {
    res <- c
  }
  colnames(res) <-gsub(" -", "\\.l", colnames(res))
  rownames(res) <-gsub("Equation ", "", rownames(res))
  return(res)
}



comp_var_coefs <- function(x) all.equal(coef_to_vars (x[[1]]), t(sapply(coef(x[[2]]), function(x) x[,"Estimate"])), check.attributes=FALSE)
comp_var_logLik <- function(x, tol=deftol) all.equal(logLik(x[[1]]), as.numeric(logLik(x[[2]])), check.attributes=FALSE, tolerance=tol)

comp_var_pred <- function(x) all.equal(predict(x[[1]]), sapply(predict(x[[2]], n.ahead=5)$fcst, function(x) x[,"fcst"]),check.attributes=FALSE)
comp_var_predOld <- function(x) all.equal(sapply(tsDyn:::predictOld.VAR(x[[1]], n.ahead=5)$fcst, function(x) x[,"fcst"]), sapply(predict(x[[2]], n.ahead=5)$fcst, function(x) x[,"fcst"]),check.attributes=FALSE)
comp_var_fevd <- function(x) all.equal(sapply(fevd(x[[1]]), head,2), sapply(fevd(x[[2]]), head,2), check.attributes=FALSE)
comp_var_IRF <- function(x, ortho = FALSE) isTRUE(all.equal(irf(x[[1]], boot=FALSE, ortho = ortho)$irf, 
                                                            irf(x[[2]], boot=FALSE, ortho = ortho)$irf, check.attributes=FALSE))
comp_var_IRF_old <- function(x, ortho = FALSE) isTRUE(all.equal(tsDyn:::irf_old(x[[1]], boot=FALSE, ortho = ortho)$irf, 
                                                                irf(x[[2]], boot=FALSE, ortho = ortho)$irf, check.attributes=FALSE))

### Compare VAR methods:
sapply(all_var_models, comp_var_coefs)
sapply(all_var_models, comp_var_logLik, tol=lowlowtol)
roundAll.Equal(sapply(all_var_models_noNoBo, comp_var_pred),7)
roundAll.Equal(sapply(all_var_models_noNoBo, comp_var_predOld),7)
sapply(all_var_models, comp_var_IRF)
sapply(all_var_models, comp_var_IRF_old)
sapply(all_var_models, comp_var_IRF, ortho = TRUE)
sapply(all_var_models, comp_var_IRF_old, ortho = TRUE)


################################
#'## Tidy
################################
comp_tidy <- function(x) {
  n_co <- x[[1]]$nparB
  
  ## need to slice because of issue: https://github.com/tidymodels/broom/issues/1174
  
  out_tsD <- tidy(x[[1]]) |> 
    as_tibble() |> 
    dplyr::mutate(term = stringr::str_replace(term, "Intercept", "const") |> 
                    stringr::str_replace("Trend", "trend") |> 
                    stringr::str_replace(" -([0-9]+)", ".l\\1")) |> 
    dplyr::rename(group=equation) |> 
    dplyr::slice(1:n_co) |> 
    dplyr::arrange(group, term)
  
  out_vars <- broom::tidy(x[[2]]) |> #as_tibble() |> 
    dplyr::slice(1:n_co) |> 
    dplyr::arrange(group, term)
  all.equal(out_vars, out_tsD)
  # waldo::compare(out_tsD,
                 # out_vars, tolerance = 1e-12)
  
}
## comp_tidy(x=all_var_models[[1]])
same_B <- sapply(all_var_models, comp_var_coefs)
sapply(all_var_models[same_B=="TRUE"], comp_tidy)

