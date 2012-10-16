

########### VAR Representation

VARrep <- function(vecm) {

  lag <- vecm$lag
  k <- vecm$k
  r <- vecm$model.specific$r
  co <- vecm$coefficients
  include <- vecm$include
  LRinclude <- vecm$model.specific$LRinclude

##obtain matrices
  betas <- vecm$model.specific$beta
  if(LRinclude!="none") betas <- betas[1:k,, drop=FALSE]

  ## Pi matrix
  Pi <-  co[, grep("ECT", colnames(co))]%*%t(betas)
  ## A_i matrix
  Amat <- matrix(NA, nrow=k, ncol=lag*(k)+k)

  if(lag>0){
    ## A_lag+1 matrix
    Amat[,(1:k)+k*lag] <- -co[,grep(paste("-", lag, sep=""), colnames(co))]
    ## A_lag+1 matrix
    if(lag>1) for(i in 1:(lag-1)) Amat[,(1:k)+k*i] <- -(co[,grep(paste("-", i, sep=""), colnames(co))] -co[,grep(paste("-", i+1, sep=""), colnames(co))])

    cumulMat <- matrix(0, k,k)
    for(i in 1:lag) cumulMat <- cumulMat + Amat[,(1:k)+k*i]
    Amat[, 1:k] <- Pi + (diag(k)- cumulMat )
  } else {
    Amat[, 1:k] <- Pi + diag(k)
  }
## Names
  varNames <- colnames(vecm$model)[1:k]
  colnames(Amat) <- paste(rep(varNames, lag+1), rep(1:(lag+1), each=k), sep=".l")

## Add deterministic terms
  if(include!="none"){
    incName <- switch(vecm$include, "const"="Intercept", trend="Trend", both="Intercept|Trend")
    incVar <- co[,grep(incName , colnames(co)),drop=FALSE]
    if(LRinclude!="none"){
      Pi_all <-  co[, grep("ECT", colnames(co))]%*%t(vecm$model.specific$beta)
      Pi_deter <- Pi_all[,"trend", drop=FALSE]
      colnames(Pi_deter) <- "Trend"
      Amat <- cbind(Pi_deter,Amat)
    }
    Amat <- cbind(incVar,Amat)
    colnames(Amat) <- gsub("Intercept", "constant", colnames(Amat))
    
  } else if(LRinclude!="none"){
    Pi_all <-  co[, grep("ECT", colnames(co))]%*%t(vecm$model.specific$beta)
    Pi_deter <- Pi_all[,switch(LRinclude, "const"="const", "trend"="trend", "both"=c("const", "trend")), drop=FALSE]
    colnames(Pi_deter) <- switch(LRinclude, "const"="constant", "trend"="Trend", "both"=c("constant", "Trend"))
    Amat <- cbind(Pi_deter,Amat)
  }
  rownames(Amat) <- gsub("Equation ","",rownames(co))
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
  rownames(co) <- gsub("Equation ", "", rownames(co))
  colnames(co) <- gsub(" -([0-9]+)","\\.l\\1", colnames(co))
  colnames(co) <- gsub("Intercept","constant", colnames(co))

## detcoeffs
  detcoeffs <- co[,grep("constant|Trend", colnames(co)), drop=FALSE]

## A
  A <- list()
  for(i in 1:lag) A[[i]] <- co[,grep(paste("\\.l", i, sep=""), colnames(co)), drop=FALSE]
  names(A) <- paste("A", 1:lag, sep="")

## Rank
  rank <- if(model=="VECM") x$model.specific$r else K

## vecm
  ecdet <- if(inherits(x,"VAR")) "none" else x$model.specific$LRinclude
  vecm<- new("ca.jo", season = NULL, dumvar=NULL, ecdet=ecdet,lag=as.integer(lag),spec="transitory")

## datamat
  if(model=="VAR"){
    datamat <- as.matrix(tail(as.data.frame(x$model),-lag))
  } else {
    newx <- lineVar(x$model[,1:K], lag=lag, include=x$include)
    datamat <- as.matrix(tail(as.data.frame(newx$model),-lag))
  }
  colnames(datamat) <- gsub(" -([0-9]+)","\\.l\\1", colnames(datamat))
  colnames(datamat) <- gsub("Intercept","constant", colnames(datamat))

## residuals
  resids <- residuals(x)
  colnames(resids) <- paste("resids of", colnames(resids))
## Return:
  result <- list(deterministic = detcoeffs, A = A, p = lag, K = K, y = as.matrix(x$model[,1:x$k]), obs = x$t, totobs = 
		  x$T, call = match.call(), vecm = vecm, datamat = datamat, resid = resids, r = rank)

  class(result) <- "vec2var"
  return(result)   

}


############################################################
#################### Methods
############################################################

predict.VAR <- function(object,...){
 predict(vec2var.tsDyn(object), ...)
}

predict.VECM <- function(object,...){
 predict(vec2var.tsDyn(object), ...)
}

irf.VAR <- function(x, impulse=NULL, response=NULL, n.ahead=10, ortho=TRUE, cumulative=FALSE, boot=TRUE, ci=0.95, runs=100, seed=NULL, ...){
 irf(vec2var.tsDyn(x), impulse=impulse, response=response, n.ahead = n.ahead, ortho=ortho, cumulative=cumulative, boot=boot, ci=ci, runs=runs, seed=seed, ...)
}

irf.VECM <- function(x, impulse=NULL, response=NULL, n.ahead=10, ortho=TRUE, cumulative=FALSE, boot=TRUE, ci=0.95, runs=100, seed=NULL, ...){
 irf(vec2var.tsDyn(x), impulse=impulse, response=response, n.ahead = n.ahead, ortho=ortho, cumulative=cumulative, boot=boot, ci=ci, runs=runs, seed=seed, ...)
}

fevd.VAR <- function(x, n.ahead=10, ...){
 fevd(vec2var.tsDyn(x),n.ahead=n.ahead, ...)
}

fevd.VECM <- function(x, n.ahead=10, ...){
 fevd(vec2var.tsDyn(x), n.ahead=n.ahead, ...)
}


predict2 <- function(vecm, n.ahead){
  lag <- vecm$lag
  k <- vecm$k
  B <- VARrep(vecm)
  original.data <- vecm$model[,1:k]
  starting <- if(is.data.frame(original.data)|is.matrix(original.data)) tail(original.data,lag+1) else if(is.ts(original.data)) apply(original.data,2,head,lag+1)
  innov <- matrix(0, nrow=n.ahead+lag+1, ncol=k)
  res <- VAR.sim(B=B, lag=lag+1, n=n.ahead+lag+1, starting=starting, innov=innov)
  colnames(res) <- colnames(original.data )
  tail(res, n.ahead)
}

predict2(vec1, n.ahead=5)


############################################################
#################### EXAMPLES, tests
############################################################

if(FALSE){


library(tsDyn)

data(zeroyld)
vec1 <- VECM(zeroyld, lag=2, estim="ML")
predict(vec1 )
fevd(vec1 )
irf(vec1 )




### Comparisons
library(vars)
data(Canada)

VECM_tsD <- VECM(Canada, lag=2, estim="ML")
VAR_tsD <- lineVar(Canada, lag=2)
VAR_tsD_tovars <-tsDyn:::vec2var.tsDyn(VAR_tsD)
VECM_tsD_tovars <-tsDyn:::vec2var.tsDyn(VECM_tsD)

VAR_vars <- VAR(Canada, p=2)
VECM_vars1 <- cajorls(ca.jo(Canada, K=3, spec="transitory"))
VECM_vars <- vec2var(ca.jo(Canada, K=3, spec="transitory"))


### Compare VECM methods:

### predict: OK!!
all.equal(predict(VECM_tsD)$fcst,predict(VECM_vars)$fcst)
all.equal(predict(VECM_tsD)$endog,predict(VECM_vars)$endog)

### fevd: OK!!
all.equal(fevd(VECM_tsD),fevd(VECM_vars))

### irf: OK!!
all.equal(irf(VECM_tsD, boot=FALSE)$irf,irf(VECM_vars, boot=FALSE)$irf)
all.equal(irf(VECM_tsD, boot=FALSE)$Lower,irf(VECM_vars, boot=FALSE)$Lower)

all.equal(irf(VECM_tsD, boot=TRUE,runs=2, seed=1234)$Lower, irf(VECM_vars, boot=TRUE,runs=2, seed=1234)$Lower)
all.equal(irf(VECM_tsD, boot=TRUE,runs=2, seed=1234)$Upper, irf(VECM_vars, boot=TRUE,runs=2, seed=1234)$Upper)


### Compare VECM methods:
predict(VAR_tsD)
all.equal(predict(VAR_tsD)$fcst,predict(VAR_vars)$fcst)
all.equal(predict(VECM_tsD)$endog,predict(VECM_vars)$endog)


### compare VARrep
data(denmark)
dat_examp <- denmark[,2:3]


toVARrep <- function(ca.jo){
  vec2<- vec2var(ca.jo)
  lags <- vec2$A[[1]]
  if(length(vec2$A)>1) for(i in 2:length(vec2$A)) lags <- cbind(lags, vec2$A[[i]])
  cbind(vec2$deterministic,lags)
}

toVARrep(ca.jo=ca.jo(dat_examp,  K=2, spec="transitory"))

VARrep(VECM(dat_examp, lag=1, include="const", estim="ML"))
toVARrep(ca.jo(dat_examp,  K=2, spec="transitory"))

all.equal(VARrep(VECM(dat_examp, lag=1, include="const", estim="ML")), toVARrep(ca.jo(dat_examp,  K=2, spec="transitory")), check.attributes=FALSE)
all.equal(VARrep(VECM(dat_examp, lag=2, include="const", estim="ML")), toVARrep(ca.jo(dat_examp,  K=3, spec="transitory")), check.attributes=FALSE)

all.equal(VARrep(VECM(dat_examp, lag=1, LRinclude="const", estim="ML")), toVARrep(ca.jo(dat_examp,  K=2, spec="transitory", ecdet="const")), check.attributes=FALSE)
all.equal(VARrep(VECM(dat_examp, lag=2, LRinclude="const", estim="ML")), toVARrep(ca.jo(dat_examp,  K=3, spec="transitory", ecdet="const")), check.attributes=FALSE)

all.equal(VARrep(VECM(dat_examp, lag=1, LRinclude="trend", estim="ML")), toVARrep(ca.jo(dat_examp,  K=2, spec="transitory", ecdet="trend")), check.attributes=FALSE)
all.equal(VARrep(VECM(dat_examp, lag=2, LRinclude="trend", estim="ML")), toVARrep(ca.jo(dat_examp,  K=3, spec="transitory", ecdet="trend")), check.attributes=FALSE)

all.equal(VARrep(VECM(Canada, lag=1, LRinclude="trend", estim="ML")), toVARrep(ca.jo(Canada,  K=2, spec="transitory", ecdet="trend")), check.attributes=FALSE)
VECM(Canada, lag=1, LRinclude="trend", estim="ML")
cajorls(ca.jo(Canada,  K=2, spec="transitory", ecdet="trend"))$rlm



#### compare slots: VECM
all.equal(VECM_tsD_tovars,VECM_vars)
all.equal(VECM_tsD_tovars$deterministic,VECM_vars$deterministic)
all.equal(VECM_tsD_tovars$A,VECM_vars$A)
all.equal(VECM_tsD_tovars$y,VECM_vars$y)
all.equal(VECM_tsD_tovars$resid,VECM_vars$resid)
attributes(VECM_vars$resid)
attributes(VECM_tsD_tovars$resid)

all.equal(VECM_tsD_tovars$datamat,VECM_vars$datamat)
all.equal(VECM_tsD_tovars$p,VECM_vars$p)
all.equal(VECM_tsD_tovars$r,VECM_vars$r)


#### compare slots: VAR


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



