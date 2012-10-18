


library(tsDyn)
library(vars)
data(Canada)


########################
######## Models
########################

##### VECM #####

### unrestricted cons
vecm_l1_co_tsD <-VECM(Canada, lag=1, include="const", estim="ML")
vecm_l3_co_tsD <-VECM(Canada, lag=3, include="const", estim="ML")

vecm_l1_co_var <- ca.jo(Canada, K=2, ecdet="none", spec="transitory")
vecm_l3_co_var <- ca.jo(Canada, K=4, ecdet="none", spec="transitory")

### restricted cons
vecm_l1_LRco_tsD <-VECM(Canada, lag=1, LRinclude="const", estim="ML")
vecm_l1_LRco_var <- ca.jo(Canada, K=2, ecdet="const", spec="transitory")

vecm_l3_LRco_tsD <-VECM(Canada, lag=3, LRinclude="const", estim="ML")
vecm_l3_LRco_var <- ca.jo(Canada, K=4, ecdet="const", spec="transitory")

### restricted trend
vecm_l1_LRtr_tsD <-VECM(Canada, lag=1, LRinclude="trend", estim="ML")
vecm_l1_LRtr_var <- ca.jo(Canada, K=2, ecdet="trend", spec="transitory")

vecm_l3_LRtr_tsD <-VECM(Canada, lag=3, LRinclude="trend", estim="ML")
vecm_l3_LRtr_var <- ca.jo(Canada, K=4, ecdet="trend", spec="transitory")

all_models <- list(
		    list(vecm_l1_co_var, vecm_l1_co_tsD), 
		    list(vecm_l3_co_var, vecm_l3_co_tsD), 
		    list(vecm_l1_LRco_var, vecm_l1_LRco_tsD),
		    list(vecm_l3_LRco_var, vecm_l3_LRco_tsD),
		    list(vecm_l1_LRtr_var, vecm_l1_LRtr_tsD),
		    list(vecm_l3_LRtr_var, vecm_l3_LRtr_tsD))

comp_teststat <- function(x) all.equal(x[[1]]@teststat, rev(rank.test(x[[2]])$res_df[,"eigen"]), check.attributes=FALSE)
comp_betas <- function(x) all.equal(cajorls(x[[1]])$beta, x[[2]]$model.specific$coint, check.attributes=FALSE)
comp_coefs <- function(x) all.equal(coefficients(cajorls(x[[1]])$rlm), t(coefficients(x[[2]])), check.attr=FALSE)
comp_LL <- function(x) all.equal(as.numeric(logLik(vec2var(x[[1]]))), logLik(x[[2]]), check.attr=FALSE)
comp_IRF <- function(x) all.equal(irf(vec2var(x[[1]]), boot=FALSE)$irf, irf(x[[2]], boot=FALSE)$irf, check.attr=FALSE)
comp_IRF_rand <- function(x) all.equal(irf(vec2var(x[[1]]), runs=2, seed=1234)$irf, irf(x[[2]], runs=2, seed=1234)$irf, check.attr=FALSE)
comp_FEVD <- function(x) all.equal(fevd(vec2var(x[[1]])), fevd(x[[2]]), check.attr=FALSE)
comp_resid <- function(x) all.equal(residuals(vec2var(x[[1]])), residuals(x[[2]]), check.attr=FALSE)
comp_fitted <- function(x) all.equal(fitted(vec2var(x[[1]])), fitted(x[[2]], level="original"), check.attr=FALSE)
comp_predict <- function(x) all.equal(predict(vec2var(x[[1]]))$fcst, predict(x[[2]])$fcst, check.attr=FALSE)


### Compare VECM methods:
print(sapply(all_models, comp_teststat ))
print(sapply(all_models, comp_betas))
print(sapply(all_models, comp_coefs))
print(sapply(all_models, comp_LL))
print(sapply(all_models, comp_IRF))
print(sapply(all_models, comp_IRF_rand))
print(sapply(all_models, comp_FEVD))
print(sapply(all_models, comp_resid))
print(sapply(all_models, comp_fitted))
print(sapply(all_models, comp_predict))

