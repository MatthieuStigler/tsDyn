TVAR.gen <- function(B, n=200, lag=1, include = c("const", "trend","none", "both"),  
                      Thresh, nthresh=1, thDelay=1,  thVar=NULL, mTh=1, 
                      starting=NULL, innov, 
                      trendStart=1,
                      show.parMat=FALSE, returnStarting=FALSE,
                      round_digits=16,
                      add.regime =FALSE){
  
  ###check correct arguments
  include<-match.arg(include)
  if(!nthresh%in%c(0,1,2))
    stop("Arg nthresh should  be either 0, 1 or 2")
  if(length(Thresh)!=nthresh)
    stop("please give as many threshold values as arg nthresh")
  if(max(thDelay)>lag)
    stop("Max of thDelay should be smaller or equal to the number of lags")
  
  ## Create some variables/parameters
  p <- lag
  k <- nrow(B) 		#Number of variables
  ninc <- switch(include, "none"=0, "const"=1, "trend"=1, "both"=2)
  npar <- k*p+ninc
  
  ### check dims B
  esp <- p*k+ninc 
  if(esp*(nthresh+1)!=ncol(B)){
    stop("Matrix B badly specified: expected ", esp*(nthresh+1), " elements ( (lag*K+ n inc)* (nthresh+1) ) but has ", ncol(B), "\n" )
  }
  
  ## starting values check dims
  if(!is.null(starting)){
    if(!all(dim(as.matrix(starting))==c(p,k)))
      stop("Bad specification of starting values. Should have nrow = lag and ncol = number of variables. But is: ", dim(as.matrix(starting)), sep="")
  }
  
  # trans var
  if (length(mTh) == 1) {
    combin <- matrix(0, ncol = 1, nrow = k)
    combin[mTh, ] <- 1
  }  else {
    combin <- matrix(mTh, ncol = 1, nrow = k)
  }
    
  
  if(nthresh%in%c(1,2))
    Thresh<-round(Thresh, round_digits)
  
  ##### put coefficients vector in right form according to arg include (arg both need no modif)
  a<-NULL
  if(include=="none") {
    for(i in 0:nthresh) a <- c(a, i*(p*k+2)+c(1,2))
  } else if(include=="const") {
    for(i in 0:nthresh) a<-c(a, i*(p*k+2)+c(2))
  } else if(include=="trend") {
    for(i in 0:nthresh) a<-c(a, i*(p*k+2)+c(1))
  }
  #if (include=="both"): correction useless
  
  Bmat <- myInsertCol(B, c=a ,0)
  nparBmat <- p*k+2
  
  
  ##############################
  ###Reconstitution boot/simul
  ##############################
  
  Yb <- matrix(0, ncol=k, nrow=n+p)
  
  # fill initial values if provided
  if(!is.null(starting)){
    Yb[1:p,] <- as.matrix(starting)
  }
  
  ## Initialise threshold variable
  if(nthresh>0){
    z2 <- vector("numeric", length=n+p)
    z2[1:p] <- Yb[1:p,]%*%combin
  }
  regime <- vector("numeric", length=n+p)
  
  ## trend, innovations
  resb <- rbind(matrix(0,nrow=p, ncol=k), 
                innov)	
  trend <- 1:(n+p) #   trend<-c(NA, 1:(T-1+add))
  trend <- c(rep(0, p), trendStart+(0:(n+p-1)))
  
  ## MAIN loop:
  if(nthresh==0){
    for(i in (p+1):(n+p)){
      Yb[i,] <- rowSums(cbind(Bmat[,1], 
                              Bmat[,2]*trend[i], 
                              Bmat[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),
                              resb[i,]))
    }
  } else if(nthresh==1){
    BDown <- Bmat[, seq_len(nparBmat)]
    BUp   <- Bmat[,-seq_len(nparBmat)]
    
    for(i in (p+1):(n+p)){
      if(round(z2[i-thDelay], round_digits)<=Thresh) {
        Yb[i,]<-rowSums(cbind(BDown[,1], 
                              BDown[,2]*trend[i], 
                              BDown[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),
                              resb[i,]))
        Yb[i,]<-round(Yb[i,], round_digits)
        regime[i] <- 1
      } else {
        Yb[i,]<-rowSums(cbind(BUp[,1], BUp[,2]*trend[i], BUp[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
        Yb[i,]<-round(Yb[i,], round_digits)
        regime[i] <- 2
      }
      z2[i]<-Yb[i,]%*%combin
    }
  }  else if(nthresh==2){
    BDown <- Bmat[,seq_len(nparBmat)]
    BMiddle <- Bmat[,seq_len(nparBmat)+nparBmat]
    BUp <- Bmat[,seq_len(nparBmat)+2*nparBmat]
    
    for(i in (p+1):(n+p)){
      if(round(z2[i-thDelay], round_digits)<=Thresh[1]) {
        Yb[i,] <- rowSums(cbind(BDown[,1], BDown[,2]*trend[i], BDown[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
        regime[i] <- 1
      } else if(round(z2[i-thDelay], round_digits)>Thresh[2]) {
        Yb[i,] <- rowSums(cbind(BUp[,1], BUp[,2]*trend[i], BUp[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
        regime[i] <- 3
      } else {
        Yb[i,] <- rowSums(cbind(BMiddle[,1],BMiddle[,2]*trend[i], BMiddle[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
        regime[i] <- 2
      }
      z2[i]<-Yb[i,]%*%combin
    }
  }
  
  ## return result
  if(show.parMat)
    print(Bmat)
  res <- round(Yb, round_digits) 
  if(add.regime) res <- cbind(res, regime)
  if(!returnStarting) res <- res[-c(1:p),, drop=FALSE] 
  return(res)
}



#'Simulation of a multivariate Threshold Autoregressive model (TVAR)
#'
#'Simulate a multivariate Threshold VAR (TVAR)
#'
#'This function offers the possibility to generate series following a TVAR.
#'
#'By giving a matrix of coefficients, on can only simulate a VAR (nthresh=0) or
#'TVAR (nthresh=1 or 2). One can have a specification with constant (default),
#'trend, both or none (see arg include). Order in parameters is include/lags
#'(VECM) and include/lags/include/lags for TVECM, hence, a matrix for a TVECM
#'with 3 regimes, a const and a 2 lags would have 2 lines and 2*(1+4) columns.
#'The innovations can be given by the user (a matrix of dim nxk, here n does
#'not include the starting values!), by default it uses a multivariate normal
#'distribution, with covariance matrix specified by varcov. The starting values
#'(of dim lags x k) can be given. The user should take care for their choice,
#'since it is not sure that the simulated values will cross the threshold even
#'once.
#'
#'The matrix \sQuote{B} has to be in the form: constant, trend, lags, then
#'repeated if many regimes. In case of uncertainty, using
#'\code{who.parMat=TRUE} will print the matrix as interpreted by the function,
#'helping the user to feed the right input.
#'
#'For the bootstrap, the function resamples data from a given TVAR model generated by
#' \code{\link{TVAR}}, returning the resampled data. 
#' A residual recursive bootstrap is used, where one uses either a simple
#' resampling, or the Wild bootstrap, either with a normal distribution (wild1) or
#' inverting the sign randomly (wild2)
#'
#'@aliases TVAR.sim
#'@param B Matrix of coefficients to simulate
#'@param Thresh The threshold value(s). Vector of length nthresh
#'@param nthresh number of threshold (see details)
#'@param n Number of observations to create when type="simul"
#'@param lag Number of lags to include in each regime
#'@param include Type of deterministic regressors to include. NOT WORKING
#'PROPERLY CURRENTLY if not const
#'@param thDelay 'time delay' for the threshold variable (as multiple of
#'embedding time delay d) PLEASE NOTE that the notation is currently different
#'to univariate models in tsDyn. The left side variable is taken at time t, and
#'not t+1 as in univariate cases.
#'@param mTh combination of variables with same lag order for the transition
#'variable. Either a single value (indicating which variable to take) or a
#'combination
#'@param starting Starting values (matrix of dimension lag x k). If not given,
#'set to zero.
#'@param innov Innovations used for simulation. Should be matrix of dim n x k.
#'By default multivariate normal. For the bootstrap case \code{TVAR.boot}, residuals are
#'resampled if argument is missing. 
#'@param varcov Variance-covariance matrix for the innovations. By default
#'identity matrix.
#'@param show.parMat Logical. Should the parameter matrix be shown? Useful to
#'understand how to give right input
#'@param \dots Further arguments passed to the underlying (un-exported)
#'\code{TVAR.gen} function
#'@return A matrix with the simulated/bootstraped series.
#'@author Matthieu Stigler
#'@seealso \code{\link{TVAR}} to estimate the TVAR.  Similar \code{\link{TVECM.sim}} and \code{\link{TVECM.boot}} for \code{\link{TVECM}}, 
#'\code{\link{VAR.sim}} and \code{\link{VAR.boot}} for VAR models estimated with \code{\link{lineVar}}. 
#'@keywords ts bootstrap
#'@export
#'@examples
#'
#'
#'## TVAR.sim: Simulation of a TVAR with 1 threshold
#'B <- rbind(c(0.11928245, 1.00880447, -0.009974585, -0.089316, 0.95425564, 0.02592617),
#'         c(0.25283578, 0.09182279,  0.914763741, -0.0530613, 0.02248586, 0.94309347))
#'colnames(B) <- paste(rep(c("Const", "Lag_1_var1", "Lag_1_var2"), 2), c("Low", "High"), sep="_")
#'sim <- TVAR.sim(B=B,nthresh=1,n=500, mTh=1, Thresh=5, starting=matrix(c(5.2, 5.5), nrow=1))
#'
#'#estimate the new serie
#'TVAR(sim, lag=1, dummyToBothRegimes=TRUE)
#'
#'
#'
TVAR.sim <- function(B, Thresh, nthresh=1, n=200, lag=1, include = c("const", "trend","none", "both"),  thDelay=1,  
                     mTh=1, starting=NULL, innov=rmnorm(n, varcov=varcov), 
                     varcov=diag(1,nrow(B)), show.parMat=FALSE,  ...){

	    TVAR.gen(B=B, Thresh=Thresh, nthresh=nthresh, n=n, lag=lag, include = include,  thDelay=thDelay,
	             mTh=mTh, starting=starting, innov=innov, show.parMat=show.parMat,  ...)
}




#'@rdname TVAR.sim
#'@param TVARobject Object of class \code{TVAR} generated by function
#'\code{\link{TVAR}}
#'@param boot.scheme The bootstrap scheme.  
#'@param seed Optional. Seed for the random resampling function.
#'@keywords ts bootstrap
#'@export
#'@examples
#'
#'
#'
#'## TVAR.boot: Bootstrap a TVAR with two threshold (three regimes)
#'data(zeroyld)
#'serie <- zeroyld
#'mod <- TVAR(data=serie,lag=1, nthresh=1)
#'TVAR.boot(mod)


TVAR.boot <- function(TVARobject, innov, seed, 
                       boot.scheme=c("resample", "wild1", "wild2", "check"), ...){
  
  if(!inherits(TVARobject, "TVAR")) stop("Please provide an object of class 'TVAR' generated by TVAR()")
  if(attr(TVARobject, "transVar")=="external" | TVARobject$model.specific$oneMatrix) {
    stop("Cannot (yet) bootstrap model with external thVar or commonInter")
  }
  boot.scheme <- match.arg(boot.scheme)
  
  t <- TVARobject$t
  k <- TVARobject$k
  lags <- TVARobject$lag

  ## retrieve starting values
  YX <- TVARobject$model
  yorig <- YX[,1:k]
  starts <- yorig[1:lags,, drop=FALSE]
  
  ## get innov, boot it
  resids <- residuals(TVARobject)
  if(!missing(seed)) set.seed(seed) 
  innov <- switch(boot.scheme, 
                  "resample"=  resids[sample(seq_len(t), replace=TRUE),], 
                  "wild1"=resids+rnorm(t), 
                  "wild2"=resids+sample(c(-1,1), size = t, replace=TRUE),
                  "check"=  resids)
  
  
  res <- TVAR.gen(B=TVARobject$coeffmat, n=t, lag=lags, include = TVARobject$include,  
                   starting=starts, innov=innov,
                   Thresh=getTh(TVARobject), thDelay=TVARobject$model.specific$thDelay,
                   nthresh=TVARobject$model.specific$nthresh,
                   mTh=TVARobject$model.specific$transCombin,
                  returnStarting=TRUE,
                   ...)
                   # show.parMat=FALSE)
  colnames(res)[1:k] <- colnames(yorig)
  res
}

# VAR.sim.old <- function(B,  n=200, lag=1, include = c("const", "trend","none", "both"),  starting=NULL,
#                      innov=rmnorm(n,  varcov=varcov), varcov=diag(1,nrow(B)), show.parMat=FALSE, round=FALSE, seed, ...){
#   TVAR.gen(B=B, nthresh=0,  n=n, lag=lag, include = include,   starting=starting, innov=innov, varcov=varcov, show.parMat=show.parMat, round=round, seed=seed, ...)
# }
# 
# VAR.boot.old <- function(VARobject, innov, seed, ...){
#   if(!inherits(VARobject, "VAR")) stop("Please provide an object of class 'VAR' generated by lineVar()")
#   if(attr(VARobject, "varsLevel")!="level") stop("Does not work with VAR in diff or ADf specification")
# 
#   TVAR.gen(TVARobject=VARobject, innov = innov, seed=seed, ...)
# }


if(FALSE){
library(tsDyn)
environment(TVAR.gen)<-environment(star)
environment(TVAR.sim)<-environment(star)
environment(TVAR.boot)<-environment(star)

##Simulation of a TVAR with 1 threshold
B<-rbind(c(0.11928245, 1.00880447, -0.009974585, -0.089316, 0.95425564, 0.02592617),c(0.25283578, 0.09182279,  0.914763741, -0.0530613, 0.02248586, 0.94309347))
sim<-TVAR.sim(B=B,nthresh=1,n=500,mTh=1, Thresh=5, starting=matrix(c(5.2, 5.5), nrow=1))

#estimate the new serie
TVAR(sim, lag=1, dummyToBothRegimes=TRUE)


##Bootstrap a TVAR 
#data(zeroyld)
serie<-zeroyld

# TVAR.sim(data=serie,nthresh=0)
# all(TVAR.sim(data=serie,nthresh=0, type="check", lag=1)==serie)

##with two threshold (three regimes)
# TVAR.sim(data=serie,nthresh=2,type="boot",mTh=1, Thresh=c(7,9))

environment(TVAR.sim)<-environment(star)

##Check case 1 (B is given) ok!
ns<-nrow(serie)

comp_tvar_sim <- function(mod, serie){
  ns <- nrow(serie)
  sim_mod <- TVAR.sim(B=coef(mod), lag=mod$lag, include=mod$include,nthresh=0, n=ns-mod$lag, innov=residuals(mod), starting=serie[1:mod$lag,,drop=FALSE])
  all.equal(sim_mod, as.matrix(serie)[-c(1:mod$lag),], check.attributes=FALSE)
}


#data(barry)
var_l1_co <-lineVar(barry, lag=1, include="const")
var_l1_tr <-lineVar(barry, lag=1, include="trend")
var_l1_bo <-lineVar(barry, lag=1, include="both")
var_l1_no <-lineVar(barry, lag=1, include="none")

var_l3_co <-lineVar(barry, lag=3, include="const")
var_l3_tr <-lineVar(barry, lag=3, include="trend")
var_l3_bo <-lineVar(barry, lag=3, include="both")
var_l3_no <-lineVar(barry, lag=3, include="none")

var_all <- list(
		var_l1_co, var_l1_tr, var_l1_bo, var_l1_no, 
		var_l3_co, var_l3_tr, var_l3_bo, var_l3_no)

lapply(var_all , comp_tvar_sim, serie=barry)
comp_tvar_sim(var_l3_tr, serie=barry)




##Check case 2 (data is given) ok!
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="const", round=FALSE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="const", lag=3, round=FALSE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="trend", round=FALSE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="trend", lag=3, round=FALSE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="both", round=FALSE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="both", lag=3, round=FALSE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="both", lag=2, round=FALSE),as.matrix(serie), check.attributes=FALSE)

all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=1, type="check",mTh=1, include="const", round=TRUE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=1, type="check",mTh=1, include="trend", round=TRUE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=1, type="check",mTh=2, include="const", round=TRUE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=1, type="check",mTh=2, include="trend", round=TRUE),as.matrix(serie), check.attributes=FALSE)

all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=2, type="check",mTh=1, include="const", round=TRUE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=2, type="check",mTh=1, include="trend", round=TRUE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=2, type="check",mTh=2, include="const", round=TRUE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=2, type="check",mTh=2, include="trend", round=TRUE),as.matrix(serie), check.attributes=FALSE)



### Check case 3: TVARobject is given
all.equal(tsDyn:::TVAR.gen(TVARobject=lineVar(serie, lag=1),type="check"),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(TVARobject=lineVar(serie, lag=1, include="trend"),type="check"),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(TVARobject=lineVar(serie, lag=1, include="both"),type="check"),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(TVARobject=lineVar(serie, lag=1, include="none"),type="check"),as.matrix(serie), check.attributes=FALSE)


all.equal(tsDyn:::TVAR.gen(TVARobject=TVAR(serie, nthresh=1, lag=1, trace=FALSE),type="check"),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(TVARobject=TVAR(serie, nthresh=1, lag=2, trace=FALSE),type="check"),as.matrix(serie), check.attributes=FALSE)

all.equal(tsDyn:::TVAR.gen(TVARobject=TVAR(serie, nthresh=2, lag=1, trace=FALSE),type="check"),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(TVARobject=TVAR(serie, nthresh=2, lag=2, trace=FALSE),type="check"),as.matrix(serie), check.attributes=FALSE)

a1 <- VAR.boot(VARobject=lineVar(serie, lag=1))
a2 <- VAR.boot(VARobject=lineVar(serie, lag=2))
a3 <- TVAR.boot(TVARobject=TVAR(serie, nthresh=1, lag=1, trace=FALSE))

}
