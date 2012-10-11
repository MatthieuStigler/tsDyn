

rank.select <- function(data, lag.max=10, r.max=ncol(data)-1, include="intercept", fitMeasure=c("SSR", "LL")) {

  fitMeasure <- match.arg(fitMeasure)

  models_list <- list()

  VARtype <- if(r.max==0) "level" else "diff"

  for(i in 1:lag.max){
    models_list[[i]] <- list()
    data_cut <- if(i==lag.max) data else data[-c(1:(lag.max-i)),]

  ## VAR: rank 0 (on diffs)
    models_list[[i]][[1]] <- try(lineVar(data_cut, lag=i, I=VARtype ), silent=TRUE)
    if(inherits(models_list[[i]][[1]], "try-error")) models_list[[i]][[1]] <- NA

  ## VECM: rank 1 to k-1
    if(r.max>0){
      for(j in 1:r.max){
	models_list[[i]][[j+1]] <- try(VECM(data_cut, lag=i, r=j, estim="ML"), silent=TRUE)
	if(inherits(models_list[[i]][[j+1]], "try-error")) models_list[[i]][[j+1]] <- FALSE
      }
      ## VAR: full rank
      models_list[[i]][[r.max+2]] <- try(lineVar(data_cut, lag=i+1, I="level"), silent=TRUE)
    }
  }

  logLik.logical <- function(x) NA


  AICs <- matrix(sapply(models_list, function(x) sapply(x,function(x) if(inherits(x, "VAR")) AIC(x,fitMeasure=fitMeasure) else NA)), ncol=lag.max)
  HCs <- matrix(sapply(models_list, function(x) sapply(x,function(model) AIC(model, k=2*log(log(model$t)),fitMeasure=fitMeasure))), ncol=lag.max)
  LLs <- matrix(sapply(models_list, function(x) sapply(x,logLik)), ncol=lag.max)
  BICs <- matrix(sapply(models_list, function(x) sapply(x,BIC, fitMeasure=fitMeasure)), ncol=lag.max)

  names_mats <- list(paste("r", 0:(r.max+min(r.max,1)), sep="="), paste("lag", 1:lag.max, sep="="))
  dimnames(AICs) <- dimnames(BICs) <- dimnames(LLs) <- dimnames(HCs)<- names_mats

## Best values:
  AIC_min <- which(AICs==min(AICs),arr.ind=TRUE)
  BIC_min <- which(BICs==min(BICs),arr.ind=TRUE)
  HC_min <- which(HCs==min(HCs),arr.ind=TRUE)

## Best rank:
  best_ranks1 <- sapply(list(AIC_min, BIC_min, HC_min), rownames)
  best_ranks <- as.numeric(gsub("r=", "", best_ranks1))
  names(best_ranks) <- c("AIC", "BIC", "HC")

## return result
  res <- list(AICs=AICs, BICs=BICs, HCs=HCs, AIC_min=AIC_min, HC_min=HC_min, BIC_min=BIC_min, LLs=LLs, best_ranks=best_ranks)
  class(res ) <- "rank.select"
  return(res)

}


lag.select <- function(data, lag.max=10, include="intercept", fitMeasure=c("SSR", "LL")) {
  rank.select(data=data, lag.max=lag.max, r.max=0, include=include, fitMeasure=fitMeasure) 
}

print.rank.select <- function(x,...){

  AIC_rank_info <- if(nrow(x$AICs)>1) paste("rank=",x$AIC_min[1,1]-1)
  BIC_rank_info <- if(nrow(x$AICs)>1) paste("rank=",x$BIC_min[1,1]-1)
  HC_rank_info <- if(nrow(x$AICs)>1) paste("rank=",x$HC_min[1,1]-1)

  cat("Best AIC:", AIC_rank_info,  " lag=", x$AIC_min[1,2],  "\n")
  cat("Best BIC:", BIC_rank_info, " lag=", x$BIC_min[1,2],  "\n")
  cat("Best HC :",  HC_rank_info, " lag=", x$HC_min[1,2],  "\n")
}

summary.rank.select <- function(x,...){

  print.rank.select(x)

  cat("\nBest number of lags:\n")
  AIC_minlag<-apply(x$AICs, 1, which.min)
  BIC_minlag<-apply(x$BICs, 1, which.min)
  HC_minlag<-apply(x$BICs, 1, which.min)

  mat <- rbind(AIC_minlag, BIC_minlag,HC_minlag)
  dimnames(mat) <- list(c("AIC", "BIC", "HC"), paste("r", if(nrow(x$BICs)==1) 0 else 0:(nrow(x$BICs)-1), sep="="))
  print(mat)
}


if(FALSE){
library(tsDyn)
library(vars)
data(Canada)

resu <- rank.select(Canada)
resu
summary(resu)


resvar_SSR <- lag.select(Canada, fitMeasure="SSR")
resvar_SSR
summary(resvar_SSR)
resvar_LL <- lag.select(Canada, fitMeasure="LL")


resvar2 <- VARselect(Canada)
resvar2$criteria[3,]
resvar_SSR$BICs/(nrow(Canada)-10)
resvar_LL$BICs/(nrow(Canada)-10)

AIC(lineVar(Canada, lag=1), fitMeasure="LL")
AIC(lineVar(Canada, lag=1), fitMeasure="SSR")


v<- lineVar(Canada, lag=1)
# AIC(v, fitMeasure="LL")
AIC(v)




#### PHILIPS DGP

### DGP r=0
n<- 200
dgp_r0_1 <- cbind(cumsum(rnorm(n)),cumsum(rnorm(n)))

### DGP r=1
alpha_1 <- matrix(c(1,0.5),ncol=1)
beta_1 <- matrix(c(-1,1),ncol=1)
Pi_1 <- alpha_1%*%t(beta_1)
# qr(Pi_1)$rank

dgp_r1 <-  VECM.sim(B=rbind(c(0.5,0,0), c(1,0,0)), beta=1, include="none",show.parMat=TRUE)
# ts.plot(dgp_r1 )

### DGP r=2
PI_r2a <- matrix(c(-0.5, 0.1, 0.2, -0.4), ncol=2, byrow=TRUE)
PI_r2b <- matrix(c(-0.5, 0.1, 0.2, -0.15), ncol=2, byrow=TRUE)
dgp_r2_a <- VAR.sim(B=PI_r2a, n=200, lag=1, include="none")
dgp_r2_b <- VAR.sim(B=PI_r2b, n=200, lag=1, include="none")

rank.select(dgp_r0_1 )
rank.select(dgp_r1 )
rank.select(dgp_r2_a )
rank.select(dgp_r2_b )



###
logLik(VECM(Canada, lag=2, r=1, estim="ML"),r=0)
logLik(lineVar(Canada, lag=2, I="diff"))

logLik(VECM(Canada, lag=2, r=1, estim="ML"), r=1)

logLik(VECM(Canada, lag=2, r=1, estim="ML"),r=2)


AIC(VECM(Canada, lag=2, r=1, estim="ML"),r=0, fitMeasure="LL")
AIC(lineVar(Canada, lag=2, I="diff"), fitMeasure="LL")


## remember: ADF(1) = diff(2) !!
deviance(linear(x=sunspots, m=2, type="level"))
deviance(linear(x=sunspots, m=1, type="ADF"))

deviance(lineVar(Canada, lag=2, I="level"))
deviance(lineVar(Canada, lag=1, I="ADF"))

}
