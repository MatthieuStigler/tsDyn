

########### VAR Representation

VARrep <- function(vecm) {

  lag <- vecm$lag
  k <- vecm$k
  r <- vecm$model.specific$r
  co <- vecm$coefficients

##obtain matrices
  ## Pi matrix
  Pi <-  co[, grep("ECT", colnames(co))]%*%t(vecm$model.specific$beta)
  ## A_i matrix
  Amat <- matrix(NA, nrow=k, ncol=lag*(k)+k)
  ## A_lag+1 matrix
  Amat[,(1:k)+k*lag] <- -co[,grep(paste("-", lag, sep=""), colnames(co))]
  ## A_lag+1 matrix
  if(lag>1) for(i in 1:(lag-1)) Amat[,(1:k)+k*i] <- -(co[,grep(paste("-", i, sep=""), colnames(co))] -co[,grep(paste("-", i+1, sep=""), colnames(co))])

  cumulMat <- matrix(0, k,k)
  for(i in 1:lag) cumulMat <- cumulMat + Amat[,(1:k)+k*i]
  Amat[, 1:k] <- Pi + (diag(k)- cumulMat )

## Names
  rownames(Amat) <- rownames(co)
  varNames <- colnames(vecm$model)[1:k]
  colnames(Amat) <- paste(rep(varNames, lag+1), rep(1:(lag+1), each=k), sep=" -")

## Add deterministic terms
  if(vecm$include!="none"){
    incName <- switch(vecm$include, "const"="Intercept", trend="Trend", both="Intercept|Trend")
    Amat <- cbind(co[,grep(incName , colnames(co)),drop=FALSE],Amat)
  }

## res
  Amat
}



############################################################
#################### vec2var.tsDyn 
############################################################



vec2var.tsDyn <- function(x){

  model <- if(inherits(x,"VECM")) "VECM" else "VAR"
  co <- coef(x)
  lag <- ifelse(model=="VECM",x$lag+1, x$lag)
  K <- x$k

## Take vec2var representation for VECMs
  if(model=="VECM") {
    co <- VARrep(x)
  }

## detcoeffs
  if(x$include=="none") warning("Not sure implemented\n")
  detcoeffs <- co[,grep("Intercept|Trend", colnames(co)), drop=FALSE]

## A
  A <- list()
  for(i in 1:lag) A[[i]] <- co[,grep(paste("-", i, sep=""), colnames(co)), drop=FALSE]
  names(A) <- paste("A", 1:lag, sep="")

## Rank
  rank <- if(model=="VECM") x$model.specific$r else K

## vecm
  setClass("vecm", representation(season="ANY", dumvar="ANY",ecdet="character", lag="integer", spec="character"))
  ecdet <- if(inherits(x,"VAR")) "none" else x$model.specific$LRinclude
  vecm<- new("vecm", season = NULL, dumvar=NULL, ecdet=ecdet,lag=as.integer(lag),spec="transitory")

## datamat
  if(model=="VAR"){
    datamat <- as.matrix(tail(as.data.frame(x$model),-lag))
  } else {
    newx <- lineVar(x$model[,1:K], lag=lag, include=x$include)
    datamat <- as.matrix(tail(as.data.frame(newx$model),-lag))
  }
## Return:
  result <- list(deterministic = detcoeffs, A = A, p = lag, K = K, y = x$model[,1:x$k], obs = x$t, totobs = 
		  x$T, call = match.call(), vecm = vecm, datamat = datamat, resid = residuals(x), r = rank)

  class(result) <- "vec2var"
  return(result)   

}


############################################################
#################### Methods
############################################################

predict.nlVar <- function(object,...){
 predict(vec2var.tsDyn(object))
}

############################################################
#################### EXAMPLES, tests
############################################################

if(FALSE){

library(vars)
library(tsDyn)

data(denmark)
dat_examp <- denmark[,2:3]

vec2var.tsDyn(VECM(dat_examp, lag=1, include="const", estim="ML"))
  vec1 <- vec2var(ca.jo(dat_examp,  K=2, spec="transitory"))$A
  cbind(vec1$A1, vec1$A2)

vec2var.tsDyn(VECM(dat_examp, lag=2, include="const", estim="ML"))
vec2 <- vec2var(ca.jo(dat_examp,  K=3, spec="transitory"))$A
  cbind(vec2$A1, vec2$A2, vec2$A3)



data(Canada)

VECM_tsD <- VECM(Canada, lag=2, estim="ML")
VAR_tsD <- lineVar(Canada, lag=2)
VECM_tsD_tovars <-vec2var.tsDyn(VECM_tsD)

VAR_vars <- VAR(Canada, p=2)
VECM_vars1 <- cajorls(ca.jo(Canada, K=3, spec="transitory"))
VECM_vars <- vec2var(ca.jo(Canada, K=3, spec="transitory"))


#### compare slots
all.equal(VECM_tsD_tovars$y,VECM_vars$y)
all.equal(VECM_tsD_tovars$deterministic,VECM_vars$deterministic, check.attributes=FALSE)
all.equal(VECM_tsD_tovars$datamat,VECM_vars$datamat, check.attributes=FALSE)
all.equal(VECM_tsD_tovars$p,VECM_vars$p, check.attributes=FALSE)

## compare coefs
coef(VECM_tsD)
t(coef(VECM_vars1$rlm))

## compare vec2var
VECM_vars$A$A3
vec2var.tsDyn(VECM_tsD)

## compare vec2var.tsDyn
VECM_tsD_tovars$A$A3
VECM_vars$A$A3

## compare residuals
head(VECM_vars$datamat)
head(VECM_tsD_tovars$datamat)

head(residuals(VECM_tsD_tovars),2)
head(residuals(VECM_vars),2)

## compare predict
VECM_tsD_tovars$r

predict(VECM_tsD_tovars, n.ahead=5)$fcst$U
predict(VECM_vars, n.ahead=5)$fcst$U

### IRF

head(irf(VECM_tsD_tovars, boot=FALSE)$irf$U,3)
head(irf(VECM_vars, boot=FALSE)$irf$U,3)

head(irf(VECM_tsD_tovars, boot=TRUE)$Upper$U,3)
head(irf(VECM_vars, boot=TRUE)$Upper$U,3)

### FEVD
head(fevd(VECM_tsD_tovars)$U,3)
head(fevd(VECM_vars)$U,3)
}



