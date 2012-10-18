

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
  if(object$include=="none"&&object$model.specific$LRinclude=="none") stop("Does not work with include='none'")
  predict(vec2var.tsDyn(object), ...)
}

predict.VECM <- function(object,...){
  if(object$include=="none"&&object$model.specific$LRinclude=="none") stop("Does not work with include='none'")
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

####### Predict 

predict2.VAR <- function(object, newdata, n.ahead=1){
  lag <- object$lag
  k <- object$k
  include <- object$include

## get coefs
  B <- coef(object)

## setup starting values (data in y), innovations (0)
  original.data <- object$model[,1:k, drop=FALSE]
  starting <-   myTail(original.data,lag)
  innov <- matrix(0, nrow=n.ahead+lag, ncol=k)

  if(!missing(newdata)) {
    if(!inherits(newdata, c("data.frame", "matrix","zoo", "ts"))) stop("Arg 'newdata' should be of class data.frame, matrix, zoo or ts")
    if(nrow(newdata)!=lag) stop("Please provide newdata with nrow=lag")
    starting <-  newdata 
  }

## use VAR sim
  res <- VAR.sim(B=B, lag=lag, n=n.ahead+lag, starting=starting, innov=innov,include=include)

## results
  colnames(res) <- colnames(original.data )
  res <- tail(res, n.ahead)

  rownames(res) <- (nrow(original.data)+1):(nrow(original.data)+n.ahead)

  return(res)
}


predict2.VECM <- function(object, n.ahead=1, newdata){
  lag <- object$lag
  k <- object$k
  include <- object$include

## get VAR rrepresentation
  B <- VARrep(object)

## setup starting values (data in y), innovations (0)
  original.data <- object$model[,1:k]
  starting <-  myTail(original.data,lag+1) 
  innov <- matrix(0, nrow=n.ahead+lag+1, ncol=k)

  if(!missing(newdata)) {
    if(!inherits(newdata, c("data.frame", "matrix","zoo", "ts"))) stop("Arg 'newdata' should be of class data.frame, matrix, zoo or ts")
    if(nrow(newdata)!=lag+1) stop("Please provide newdata with nrow=lag+1 (note lag=p in VECM representation corresponds to p+1 in VAR rep)")
    starting <-  newdata 
  }

## use VAR sim
  res <- VAR.sim(B=B, lag=lag+1, n=n.ahead+lag+1, starting=starting, innov=innov, include=include)

## results
  colnames(res) <- colnames(original.data )
  res <- tail(res, n.ahead)
  rownames(res) <- (nrow(original.data)+1):(nrow(original.data)+n.ahead)

  return(res)
}


myHead <- function(x, n=6){

  if(inherits(x, "ts")){
    res <-  apply(x,2,head,n)
    if(n==1) res <- matrix(res, nrow=1)
  } else {
    res <-  head(x,n) 
  }

  res
}

myTail <- function(x, n=6){

  if(inherits(x, "ts")){
    res <-  apply(x,2,tail,n)
    if(n==1) res <- matrix(res, nrow=1)
  } else {
    res <-  tail(x,n) 
  }

  res
}

############################################################
#################### EXAMPLES, tests
############################################################

if(FALSE){


library(tsDyn)

data(zeroyld)
vec1 <- VECM(zeroyld, lag=2, estim="ML")
predict(vec1 )
predict2.VECM(vec1, n.ahead=5)
fevd(vec1 )
irf(vec1, runs=10 )




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

### predict


### compare prediction and actual
Var_1 <- lineVar(Canada, lag=1)
n <- nrow(Canada)
all.equal(predict2.VAR(Var_1),predict2.VAR(Var_1, newdata=Canada[n,,drop=FALSE]))
all.equal(tail(fitted(Var_1),1),predict2.VAR(Var_1, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)

Var_2 <- lineVar(Canada, lag=2)
all.equal(predict2.VAR(Var_2),predict2.VAR(Var_2, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(tail(fitted(Var_2),1),predict2.VAR(Var_2, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)

Var_1_t <- lineVar(Canada, lag=1, include="trend")
all.equal(predict2.VAR(Var_1_t),predict2.VAR(Var_1_t, newdata=Canada[n,,drop=FALSE]))
all.equal(tail(fitted(Var_1_t),1),predict2.VAR(Var_1_t, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)

Var_1_no <- lineVar(Canada, lag=1, include="none")
all.equal(predict2.VAR(Var_1_no),predict2.VAR(Var_1_no, newdata=Canada[n,,drop=FALSE]))
all.equal(tail(fitted(Var_1_no),1),predict2.VAR(Var_1_no, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)

Var_1_bo <- lineVar(Canada, lag=1, include="both")
all.equal(predict2.VAR(Var_1_bo),predict2.VAR(Var_1_bo, newdata=Canada[n,,drop=FALSE]))
all.equal(tail(fitted(Var_1_bo),1),predict2.VAR(Var_1_bo, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)


### VECM
Vecm_1_co <- VECM(Canada, lag=1, include="const")
all.equal(predict2.VECM(Vecm_1_co),predict2.VECM(Vecm_1_co, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(predict2.VECM(Vecm_1_co), varpToDf( predict(Vecm_1_co,n.ahead=1)), check.attributes=FALSE)
all.equal(tail(fitted(Vecm_1_co, level="original"),1),predict2.VECM(Vecm_1_co, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)


Vecm_2_co <- VECM(Canada, lag=2, include="const")
all.equal(predict2.VECM(Vecm_2_co),predict2.VECM(Vecm_2_co, newdata=Canada[c(n-2,n-1,n),,drop=FALSE]))
all.equal(predict2.VECM(Vecm_2_co), varpToDf( predict(Vecm_2_co,n.ahead=1)), check.attributes=FALSE)
all.equal(tail(fitted(Vecm_2_co, level="original"),1),predict2.VECM(Vecm_2_co, n.ahead=1, newdata=Canada[c(n-3,n-2,n-1),,drop=FALSE]), check.attributes=FALSE)

Vecm_1_tr <- VECM(Canada, lag=1, include="trend")
all.equal(predict2.VECM(Vecm_1_tr),predict2.VECM(Vecm_1_tr, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(predict2.VECM(Vecm_1_tr), varpToDf( predict(Vecm_1_tr,n.ahead=1)), check.attributes=FALSE)
all.equal(tail(fitted(Vecm_1_tr, level="original"),1),predict2.VECM(Vecm_1_tr, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)

Vecm_1_non <- VECM(Canada, lag=1, include="none")
all.equal(predict2.VECM(Vecm_1_non),predict2.VECM(Vecm_1_non, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(tail(fitted(Vecm_1_non, level="original"),1),predict2.VECM(Vecm_1_non, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)

Vecm_1_LRco <- VECM(Canada, lag=1, LRinclude="const")
all.equal(predict2.VECM(Vecm_1_LRco),predict2.VECM(Vecm_1_LRco, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(tail(fitted(Vecm_1_LRco, level="original"),1),predict2.VECM(Vecm_1_LRco, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)



Vecm_1_LRc <- VECM(Canada, lag=1, LRinclude="const")
all.equal(predict2.VECM(Vecm_1_LRc),predict2.VECM(Vecm_1_LRc, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(tail(fitted(Vecm_1_LRc),1),predict2.VECM(Vecm_1_LRc, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)

Vecm_1_LRt <- lineVar(Canada, lag=1, LRinclude="trend")
all.equal(predict2.VECM(Vecm_1_LRt),predict2.VECM(Vecm_1_LRt, newdata=Canada[n,,drop=FALSE]))
all.equal(tail(fitted(Vecm_1_LRt),1),predict2.VECM(Vecm_1_LRt, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)

Vecm_1_LRt_noc <- lineVar(Canada, lag=1, LRinclude="trend", include="none")
all.equal(predict2.VECM(Vecm_1_LRt_noc),predict2.VECM(Vecm_1_LRt_noc, newdata=Canada[n,,drop=FALSE]))
all.equal(tail(fitted(Vecm_1_LRt_noc),1),predict2.VECM(Vecm_1_LRt_noc, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)

Vecm_1_LRbo<- lineVar(Canada, lag=1, LRinclude="both")
all.equal(predict2.VECM(Vecm_1_LRbo),predict2.VECM(Vecm_1_LRbo, newdata=Canada[n,,drop=FALSE]))
all.equal(tail(fitted(Vecm_1_LRbo),1),predict2.VECM(Vecm_1_LRbo, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)

predict2.VAR(lineVar(Canada, lag=1, include="none"))
predict2.VECM(VECM(Canada, lag=1))
predict2.VECM(VECM(Canada, lag=1), newdata=Canada[(nrow(Canada)-1):nrow(Canada),,drop=FALSE])


varpToDf <- function(x) matrix(sapply(x$fcst, function(x) x[,"fcst"]), nrow=nrow(x$fcst[[1]]))

all.equal(sapply(predict(VECM_tsD, n.ahead=5)$fcst, function(x) x[,"fcst"]), predict2(VECM_tsD, n.ahead=5), check.attributes=FALSE, tol=1e-07)
all.equal(sapply(predict(VECM_tsD, n.ahead=10)$fcst, function(x) x[,"fcst"]), predict2(VECM_tsD, n.ahead=10), check.attributes=FALSE, tol=1e-07)


### VAR
all.equal(predict2.VAR(lineVar(Canada, lag=2), n.ahead=3), sapply(predict(VAR(Canada, p=2), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)
all.equal(predict2.VAR(lineVar(Canada, lag=4), n.ahead=3), sapply(predict(VAR(Canada, p=4), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)


all.equal(predict2.VAR(lineVar(Canada, lag=1, include="none"), n.ahead=3), sapply(predict(VAR(Canada, p=1, type="none"), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)
all.equal(predict2.VAR(lineVar(Canada, lag=2, include="none"), n.ahead=3), sapply(predict(VAR(Canada, p=2, type="none"), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)

all.equal(predict2.VAR(lineVar(Canada, lag=1, include="trend"), n.ahead=3), sapply(predict(VAR(Canada, p=1, type="trend"), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)
all.equal(predict2.VAR(lineVar(Canada, lag=2, include="trend"), n.ahead=3), sapply(predict(VAR(Canada, p=2, type="trend"), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)

all.equal(predict2.VAR(lineVar(Canada, lag=1, include="both"), n.ahead=3), sapply(predict(VAR(Canada, p=1, type="both"), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)
all.equal(predict2.VAR(lineVar(Canada, lag=2, include="both"), n.ahead=3), sapply(predict(VAR(Canada, p=2, type="both"), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)

### VECM
all.equal(sapply(predict(VECM_tsD, n.ahead=5)$fcst, function(x) x[,"fcst"]), predict2(VECM_tsD, n.ahead=5), check.attributes=FALSE)
all.equal(sapply(predict(VECM_tsD, n.ahead=10)$fcst, function(x) x[,"fcst"]), predict2(VECM_tsD, n.ahead=10), check.attributes=FALSE)


ve2_tsD <- VECM(Canada, lag=2, estim="ML")
ve2_var <- vec2var(ca.jo(Canada, K=3, spec="transitory"))
all.equal(sapply(predict(ve2_tsD, n.ahead=5)$fcst, function(x) x[,"fcst"]), predict2(ve2_tsD, n.ahead=5), check.attributes=FALSE)
all.equal(sapply(predict(ve2_var, n.ahead=10)$fcst, function(x) x[,"fcst"]), predict2(ve2_tsD, n.ahead=10), check.attributes=FALSE)

}



