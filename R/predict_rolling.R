
predict_rolling  <- function (object, ...)  
  UseMethod("predict_rolling")

predict_rolling.default  <- function (object, ...)  NULL



predict_rolling_1step.nlVar <- function(object, nroll=10, n.ahead=1, refit.every, newdata, ...){

## Checks
  if(!missing(refit.every)&&refit.every>nroll) stop("arg 'refit.every' should be smaller or equal to arg 'nroll'")

## infos on model
  model <- attr(object, "model")
  k <- object$k
  origSerie <- object$model[,1:k]
  lag <- object$lag
  include <- object$include
  T <- object$T

## Create Refit (modfit) function:
  if(model=="VAR"){
    level <- attr(object, "varsLevel")
    modFit <- function(dat) lineVar(dat, lag=lag, include=include, I=level)
    add <- if(level=="level") 0 else 1
  } else {
    r <- object$model.specific$r
    estim <- object$model.specific$estim
    if(estim =="OLS") estim <- "2OLS"
    LRinclude <- object$model.specific$LRinclude
    modFit <- function(dat) VECM(dat,   lag=lag, include=include, estim=estim, r=r, LRinclude=LRinclude)
    add <- 1
  }


## Set refit.every
  everys <-  if(!missing(refit.every)) seq(refit.every, by=refit.every, to=nroll) else 0

## Fit initial model
  if(missing(newdata)){
    subSerie <- myHead(origSerie, -nroll)
    outSerie <- myTail(origSerie, nroll+lag+add)
    mod <- modFit(subSerie)
  } else {
    subSerie <- origSerie
    outSerie <- newdata
    if(!isTRUE(all.equal(myHead(outSerie,lag+add),myTail(subSerie,lag+add), check.attributes=FALSE))) {
      print(myHead(outSerie,lag))
      print(myTail(subSerie,lag))
      warning("'newdata' should contain as first values the last values taken for estimation, as these will be the basis for first forecast")
      outSerie <- rbind(myHead(subSerie,lag), outSerie)
    }
    if(nrow(newdata)!=nroll+lag-1+add) {
      warning("nroll adjusted to dimension of newdata")
      nroll <- nrow(newdata)
    }
    mod <- object
  }

## Refit model on smaller sample:
  R <- matrix(0, ncol=k, nrow=nroll)
  colnames(R) <- colnames(origSerie)

  for(i in 1:nroll){

  ## model
    if(i%in%everys){
      subSerie <- myHead(origSerie, -nroll+i-1)
      mod <- modFit(subSerie)
      out <- outSerie[i:(i+lag-1+add),,drop=FALSE]
    } else {
      out <- outSerie[(i):(i+lag-1+add),,drop=FALSE]
    }

  ## pred
    R[i,] <- predict(mod, n.ahead=n.ahead, newdata=out)[n.ahead,]

    lastPos <- T-nroll-n.ahead+i
    lags <- c(0:max(0,lag-1+add))
print(lastPos-lags)
    dat <- origSerie[lastPos-lags,,drop=FALSE] # old: #     dat <- myTail(origSerie[1:(T-nroll-n.ahead+1),], lag-1+add)
    R[i,] <- predict(mod, n.ahead=n.ahead, newdata=dat)[n.ahead,]
  }


## Return
  res <- list(pred=as.data.frame(R), true=outSerie)
  return(res)

}



predict_rolling.nlVar<- function(object, nroll=10, n.ahead=1:2, refit.every, newdata, ...){

  morgAr <- list(object=object, nroll=nroll)
  if(!missing(refit.every)) morgAr$refit.every <- refit.every
  if(!missing(newdata)) morgAr$newdata <- newdata

  if(length(n.ahead)==1){
    if(missing(newdata)){
      res_predroll <- predict_rolling_1step.nlVar(object, n.ahead=n.ahead, nroll=nroll, refit.every=refit.every,...)
    } else {
      res_predroll <- predict_rolling_1step.nlVar(object, n.ahead=n.ahead, nroll=nroll, newdata=newdata, refit.every=refit.every,...)
    }
    res <- res_predroll$pred
    newdata <- res_predroll$true
  } else {
    res_map <- mapply(predict_rolling_1step.nlVar, n.ahead=n.ahead, MoreArgs=morgAr,SIMPLIFY = FALSE)
    res_li <- lapply(res_map, function(x) x$pred)
    if(missing(newdata)) newdata <- res_map[[1]]$true ## VECM case
    res <- as.data.frame(simplify2df(res_li))

    res$n.ahead <- rep(n.ahead, each=nrow(res)/length(n.ahead))
  }

## return result
  res <- list(pred=res, true=data.frame(newdata))
  return(res)
}


predict_rolling.nlar <- function(object, n.ahead=1, newdata, ...){

  if(missing(newdata)) stop("Providing newdata required for objects nlar")
  if(length(newdata) > length(object$str$x)) {
    stop("newdta should not have length bigger than sample used to estimate 'object'. Be careful not to provide first sub-sample in newdata!")
  }
#   n.aheads <- length(n.ahead)
 
## Construct data
  estim_samp <- object$str$x
  n_estim_samp <- length(estim_samp)
  nroll <- length(newdata)
  full_samp <- c(estim_samp, newdata)

  pred <- vector("numeric",length=nroll*length(n.ahead))
  for(j in 1:length(n.ahead)){
    n.ahead_i <- n.ahead[j]
    for(i in 1:nroll){
      pred[i+(j-1)*nroll] <- predict(object, n.ahead=n.ahead_i, newdat=full_samp[1:(n_estim_samp+i-n.ahead_i)])[n.ahead_i]
    }
  }
  if(length(n.ahead)>1){
    pred <- data.frame(pred=pred,n.ahead =rep(n.ahead, each=nroll))
  } else {
    pred <- as.data.frame(pred)
  }
  colnames(pred)[1] <- object$str$series

  ## Return object
  res <- list(pred=pred, true=as.data.frame(newdata))
  return(res)

}







simplify2df <- function(x) {

  out <- x[[1]]
  for(i in 2:length(x)){
    out <- rbind(out, x[[i]])
  }
  out
}


predict_rolling_fcstpkg <- function(object, n.ahead=1, newdata, model, check=FALSE, ...){

  mod_cl <- deparse(substitute(model))
  if(missing(newdata)) stop("Providing newdata required for objects ",  mod_cl, "!")
  if(length(newdata) > length(object$x)) stop("newdta should not have leength bigger than sample used to estimate 'object'. Be careful not to provide first sub-sample in newdata!")
 

  estim_samp <- object$x
  n_estim_samp <- length(estim_samp)
  nroll <- length(newdata)
  full_samp <- c(estim_samp, newdata)

## Get in sample forecasts
  slotMod<- switch(mod_cl, "Arima"="pred", "ets"="mean")
  pred <- vector("numeric",length=nroll*length(n.ahead))

  for(j in 1:length(n.ahead)){
    for(i in 1:nroll){
      mod <- model(full_samp[1:(n_estim_samp+i-n.ahead[j])], model=object)
      pred[i+(j-1)*nroll] <- forecast(mod, n.ahead=n.ahead[j])$mean[n.ahead[j]] 
    }
  }

## format pred: add eventually n.ahead column
  if(length(n.ahead)>1){
    pred <- data.frame(pred=pred,n.ahead =rep(n.ahead, each=nroll))
  } else {
    pred <- as.data.frame(pred)
  }
  colnames(pred)[1] <- if(mod_cl=="Arima") object$series else deparse(object$call$y)

## Return object
  res <- list(pred=pred, true= as.data.frame(newdata))#, model=mod_cl)
#   class(res) <- "fcstpkg"
  return(res)
}


predict_rolling.Arima <- function(object, n.ahead=1, newdata,  ...){
  predict_rolling_fcstpkg(object=object,  n.ahead=n.ahead, newdata=newdata, model=Arima,check=TRUE,  ...)
}

predict_rolling.ets <- function(object,  n.ahead=1, newdata,  ...){
  predict_rolling_fcstpkg(object=object, n.ahead=n.ahead, newdata=newdata, model=ets, check=FALSE, ...)
}


##############################################################
##################### TESTS
##############################################################
if(FALSE){
library(tsDyn)
data(barry)
n_ca<- nrow(barry)

 environment(predict_rolling_1step.nlVar) <- environment(star)
# environment(predict_rolling.nlVar) <- environment(star)
# environment(predict_rolling.default) <- environment(star)
# environment(predict_rolling) <- environment(star)

#### No refit lag=1
as.M <- function(x) as.matrix(x)
mod_var_1_full <- lineVar(barry, lag=1)
mod_var_1_sub <- lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1)


pred_roll_1<-predict_rolling_1step.nlVar(object=mod_var_1_full, nroll=10, n.ahead=1)
pred_roll_2 <-predict_rolling(object=mod_var_1_full, nroll=10, n.ahead=2)
pred_roll_12<-predict_rolling(object=mod_var_1_full, nroll=10, n.ahead=1:2)
all.equal(rbind(pred_roll_1$pred, pred_roll_2$pred), pred_roll_12$pred[,-4], check=FALSE) ## internal consistency

pred_0_12_nd <- predict(mod_var_1_sub, n.ahead=2, newdata=barry[n_ca-11,,drop=FALSE])
pred_1_12_nd <- predict(mod_var_1_sub, n.ahead=2, newdata=barry[n_ca-10,,drop=FALSE])
pred_2_12_nd <- predict(mod_var_1_sub, n.ahead=2, newdata=barry[n_ca-9,,drop=FALSE])
pred_1_12 <- predict(mod_var_1_sub, n.ahead=2)
all.equal(pred_1_12_nd, pred_1_12) ## minor: consistency in predict with/withotut newdata=dataset


pred_nd <- rbind(pred_0_12_nd, pred_1_12_nd, pred_2_12_nd)
all.equal(pred_nd[c(3,5),], as.M(pred_roll_12$pred[1:2,1:3]), check=FALSE) ## check 1-ahead
all.equal(pred_nd[c(2,4),], as.M(pred_roll_12$pred[11:12,1:3]), check=FALSE) ## check 2 ahead


#### No refit lag=3
mod_var_l3_full <- lineVar(barry, lag=3)
mod_var_l3_sub <- lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3)

pred_roll_1<-predict_rolling_1step.nlVar(object=mod_var_l3_full, nroll=10, n.ahead=1)
pred_roll_2 <-predict_rolling(object=mod_var_l3_full, nroll=10, n.ahead=2)
pred_roll_12<-predict_rolling(object=mod_var_l3_full, nroll=10, n.ahead=1:2)
all.equal(rbind(pred_roll_1$pred, pred_roll_2$pred), pred_roll_12$pred[,-4], check=FALSE) ## internal consistency

pred_0_12_nd <- predict(mod_var_l3_sub, n.ahead=2, newdata=barry[n_ca-c(13:11),,drop=FALSE])
pred_1_12_nd <- predict(mod_var_l3_sub, n.ahead=2, newdata=barry[n_ca-(12:10),,drop=FALSE])
pred_2_12_nd <- predict(mod_var_l3_sub, n.ahead=2, newdata=barry[n_ca-c(9:11),,drop=FALSE])
pred_1_12 <- predict(mod_var_l3_sub, n.ahead=2)
all.equal(pred_1_12_nd, pred_1_12) ## minor: consistency in predict with/withotut newdata=dataset


pred_nd <- rbind(pred_0_12_nd, pred_1_12_nd, pred_2_12_nd)
all.equal(pred_nd[c(3,5),], as.M(pred_roll_12$pred[1:2,1:3]), check=FALSE) ## check 1-ahead
all.equal(pred_nd[c(2,4),], as.M(pred_roll_12$pred[11:12,1:3]), check=FALSE) ## check 2 ahead


pred_roll<-predict_rolling(object=lineVar(barry, lag=3), nroll=10, n.ahead=1)

pred1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-2):n_ca)-10,,drop=FALSE])
pred2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-2):n_ca)-9,,drop=FALSE])
all.equal(rbind(pred1, pred2), head(pred_roll$pred,2), check=FALSE)

all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-2):n_ca)-1,,drop=FALSE]), tail(pred_roll$pred,1), check=FALSE)



### Refit lag=1
pred_ref<-predict_rolling(object=lineVar(barry, lag=1), nroll=10, n.ahead=1, refit.every=1)

pred_ref_1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1)
pred_ref_2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-9), lag=1), n.ahead=1)
all.equal(rbind(pred_ref_1, pred_ref_2), head(pred_ref$pred,2), check=FALSE)
all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-1), lag=1), n.ahead=1), tail(pred_ref$pred,1), check=FALSE)

### Refit lag=3
pred_ref<-predict_rolling(object=lineVar(barry, lag=3), nroll=10, n.ahead=1, refit.every=1)

pred_ref_1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1)
pred_ref_2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-9), lag=3), n.ahead=1)
all.equal(rbind(pred_ref_1, pred_ref_2), head(pred_ref$pred,2), check=FALSE)
all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-1), lag=3), n.ahead=1), tail(pred_ref$pred,1), check=FALSE)


#### No refit: VAR diff,  lag=1
pred_roll_dif<-predict_rolling(object=lineVar(barry, lag=1, I="diff"), nroll=10, n.ahead=1)

pred1_dif <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1, I="diff"), n.ahead=1, newdata=barry[(n_ca-1):n_ca-10,,drop=FALSE])
pred2_dif <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1, I="diff"), n.ahead=1, newdata=barry[(n_ca-1):n_ca-9,,drop=FALSE])
all.equal(rbind(pred1_dif, pred2_dif), head(pred_roll_dif$pred,2), check=FALSE)

all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1, I="diff"), n.ahead=1, newdata=barry[(n_ca-1):n_ca-1,,drop=FALSE]), tail(pred_roll_dif$pred,1), check=FALSE)


#### VECM No refit lag=1
pred_vec_roll_l1 <- predict_rolling(object=VECM(barry, lag=1), nroll=10, n.ahead=1)

pred_vec_l1_1 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[((n_ca-1):n_ca)-10,,drop=FALSE])
pred_vec_l1_2 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[((n_ca-1):n_ca)-9,,drop=FALSE])
all.equal(rbind(pred_vec_l1_1, pred_vec_l1_2), head(pred_vec_roll_l1$pred,2), check=FALSE)

all.equal(predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[((n_ca-1):n_ca)-1,,drop=FALSE]), tail(pred_vec_roll_l1$pred,1), check=FALSE)



#### VECM No refit lag=3
pred_vec_roll <- predict_rolling(object=VECM(barry, lag=3), nroll=10, n.ahead=1)

pred_vec1 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-3):n_ca)-10,,drop=FALSE])
pred_vec2 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-3):n_ca)-9,,drop=FALSE])
all.equal(rbind(pred_vec1, pred_vec2), head(pred_vec_roll$pred,2), check=FALSE)

all.equal(predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-3):n_ca)-1,,drop=FALSE]), tail(pred_vec_roll$pred,1), check=FALSE)



### VECM Refit lag=1
pred_vec_ref<-predict_rolling(object=VECM(barry, lag=1), nroll=10, n.ahead=1, refit.every=1)

pred_vec_ref_1 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1)
pred_vec_ref_2 <- predict(VECM(tsDyn:::myHead(barry,n_ca-9), lag=1), n.ahead=1)
all.equal(rbind(pred_vec_ref_1, pred_vec_ref_2), head(pred_vec_ref$pred,2), check=FALSE)
all.equal(predict(VECM(tsDyn:::myHead(barry,n_ca-1), lag=1), n.ahead=1), tail(pred_vec_ref$pred,1), check=FALSE)

############################
########### NLAR
############################

#### linear
mod_ar <- linear(lynx[1:100], m=1)

pr_ar_1 <- predict_rolling(mod_ar, newdata=log(lynx[101:114]), n.ahead=1)$pred
pr_ar_2 <- predict_rolling(mod_ar, newdata=log(lynx[101:114]), n.ahead=2)$pred

pr_ar_12 <- predict_rolling(mod_ar, newdata=log(lynx[101:114]), n.ahead=1:2)$pred
all.equal(c(pr_ar_1[,1],pr_ar_2[,1]), pr_ar_12[,1])
all.equal(c(pr_ar_1[1,1],pr_ar_2[2,1]), predict(mod_ar, n.ahead=3)[1:2])



#### setar
mod_set <- setar(lynx[1:100], m=1)

pr_set_1 <- predict_rolling(mod_set, newdata=log(lynx[101:114]), n.ahead=1)$pred
pr_set_2 <- predict_rolling(mod_set, newdata=log(lynx[101:114]), n.ahead=2)$pred

pr_set_12 <- predict_rolling(mod_set, newdata=log(lynx[101:114]), n.ahead=1:2)$pred
all.equal(c(pr_set_1[,1],pr_set_2[,1]), pr_set_12[,1])
all.equal(c(pr_set_1[1,1],pr_set_2[2,1]), predict(mod_set, n.ahead=3)[1:2])


############################
########### FORECATS
############################
library(forecast)
mod_arauto <- auto.arima(log(lynx[1:100]))
mod_ets <- ets(log(lynx[1:100]))
mod_arim <- Arima(log(lynx[1:100]), order=c(1,0,0))




## ARIMA
pr_fct_1 <- predict_rolling(mod_arim, newdata=log(lynx[101:114]), n.ahead=1)$pred
pr_fct_2 <- predict_rolling(mod_arim, newdata=log(lynx[101:114]), n.ahead=2)$pred

pr_fct_12 <- predict_rolling(mod_arim, newdata=log(lynx[101:114]), n.ahead=1:2)$pred
all.equal(c(pr_fct_1[,1],pr_fct_2[,1]), pr_fct_12[,1])
all.equal(c(pr_fct_1[1,1],pr_fct_2[2,1]), forecast(mod_arim, h=3)$mean[1:2])

## auto.ARIMA
pr_fct_at_1 <- predict_rolling(mod_arauto, newdata=log(lynx[101:114]), n.ahead=1)$pred
pr_fct_at_2 <- predict_rolling(mod_arauto, newdata=log(lynx[101:114]), n.ahead=2)$pred

pr_fct_at_12 <- predict_rolling(mod_arauto, newdata=log(lynx[101:114]), n.ahead=1:2)$pred
all.equal(c(pr_fct_at_1[,1],pr_fct_at_2[,1]), pr_fct_at_12[,1])
all.equal(c(pr_fct_at_1[1,1],pr_fct_at_2[2,1]), forecast(mod_arauto, h=3)$mean[1:2])

## ETS
pr_fct_ets_1 <- predict_rolling(mod_ets, newdata=log(lynx[101:114]), n.ahead=1)$pred
pr_fct_ets_2 <- predict_rolling(mod_ets, newdata=log(lynx[101:114]), n.ahead=2)$pred

pr_fct_ets_12 <- predict_rolling(mod_ets, newdata=log(lynx[101:114]), n.ahead=1:2)$pred
all.equal(c(pr_fct_ets_1[,1],pr_fct_ets_2[,1]), pr_fct_ets_12[,1])
all.equal(c(pr_fct_ets_1[1,1],pr_fct_ets_2[2,1]), forecast(mod_ets, h=3)$mean[1:2])


}
