
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

## Refit function:
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

  }


## Return
  res <- list(pred=R, true=outSerie)
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
# browser()
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




##############################################################
##################### TESTS
##############################################################
if(FALSE){
library(tsDyn)
data(barry)
n_ca<- nrow(barry)

environment(predict_rolling_1step.nlVar) <- environment(star)
environment(predict_rolling.nlVar) <- environment(star)
environment(predict_rolling.default) <- environment(star)
environment(predict_rolling) <- environment(star)

#### No refit lag=1
pred_roll<-predict_rolling(object=lineVar(barry, lag=1), nroll=10, n.ahead=1)
pred_roll_vec_1<-predict_rolling(object=lineVar(barry, lag=1), nroll=10, n.ahead=1)
pred_roll_vec_12<-predict_rolling.default(object=lineVar(barry, lag=1), nroll=10, n.ahead=1:2)

pred1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[n_ca-10,,drop=FALSE])
pred2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[n_ca-9,,drop=FALSE])
all.equal(rbind(pred1, pred2), head(pred_roll$pred,2), check=FALSE)

all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[n_ca-1,,drop=FALSE]), tail(pred_roll$pred,1), check=FALSE)

#### No refit lag=3
pred_roll<-predict_rolling.nlVar(object=lineVar(barry, lag=3), nroll=10, n.ahead=1)

pred1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-2):n_ca)-10,,drop=FALSE])
pred2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-2):n_ca)-9,,drop=FALSE])
all.equal(rbind(pred1, pred2), head(pred_roll$pred,2), check=FALSE)

all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-2):n_ca)-1,,drop=FALSE]), tail(pred_roll$pred,1), check=FALSE)



### Refit lag=1
pred_ref<-predict_rolling.nlVar(object=lineVar(barry, lag=1), nroll=10, n.ahead=1, refit.every=1)

pred_ref_1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1)
pred_ref_2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-9), lag=1), n.ahead=1)
all.equal(rbind(pred_ref_1, pred_ref_2), head(pred_ref$pred,2), check=FALSE)
all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-1), lag=1), n.ahead=1), tail(pred_ref$pred,1), check=FALSE)

### Refit lag=3
pred_ref<-predict_rolling.nlVar(object=lineVar(barry, lag=3), nroll=10, n.ahead=1, refit.every=1)

pred_ref_1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1)
pred_ref_2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-9), lag=3), n.ahead=1)
all.equal(rbind(pred_ref_1, pred_ref_2), head(pred_ref$pred,2), check=FALSE)
all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-1), lag=3), n.ahead=1), tail(pred_ref$pred,1), check=FALSE)


#### No refit: VAR diff,  lag=1
pred_roll_dif<-predict_rolling.nlVar(object=lineVar(barry, lag=1, I="diff"), nroll=10, n.ahead=1)

pred1_dif <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1, I="diff"), n.ahead=1, newdata=barry[(n_ca-1):n_ca-10,,drop=FALSE])
pred2_dif <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1, I="diff"), n.ahead=1, newdata=barry[(n_ca-1):n_ca-9,,drop=FALSE])
all.equal(rbind(pred1_dif, pred2_dif), head(pred_roll_dif$pred,2), check=FALSE)

all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1, I="diff"), n.ahead=1, newdata=barry[(n_ca-1):n_ca-1,,drop=FALSE]), tail(pred_roll_dif$pred,1), check=FALSE)


#### VECM No refit lag=1
pred_vec_roll_l1 <- predict_rolling.nlVar(object=VECM(barry, lag=1), nroll=10, n.ahead=1)

pred_vec_l1_1 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[((n_ca-1):n_ca)-10,,drop=FALSE])
pred_vec_l1_2 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[((n_ca-1):n_ca)-9,,drop=FALSE])
all.equal(rbind(pred_vec_l1_1, pred_vec_l1_2), head(pred_vec_roll_l1$pred,2), check=FALSE)

all.equal(predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[((n_ca-1):n_ca)-1,,drop=FALSE]), tail(pred_vec_roll_l1$pred,1), check=FALSE)



#### VECM No refit lag=3
pred_vec_roll <- predict_rolling.nlVar(object=VECM(barry, lag=3), nroll=10, n.ahead=1)

pred_vec1 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-3):n_ca)-10,,drop=FALSE])
pred_vec2 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-3):n_ca)-9,,drop=FALSE])
all.equal(rbind(pred_vec1, pred_vec2), head(pred_vec_roll$pred,2), check=FALSE)

all.equal(predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-3):n_ca)-1,,drop=FALSE]), tail(pred_vec_roll$pred,1), check=FALSE)



### VECM Refit lag=1
pred_vec_ref<-predict_rolling.nlVar(object=VECM(barry, lag=1), nroll=10, n.ahead=1, refit.every=1)

pred_vec_ref_1 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1)
pred_vec_ref_2 <- predict(VECM(tsDyn:::myHead(barry,n_ca-9), lag=1), n.ahead=1)
all.equal(rbind(pred_vec_ref_1, pred_vec_ref_2), head(pred_vec_ref$pred,2), check=FALSE)
all.equal(predict(VECM(tsDyn:::myHead(barry,n_ca-1), lag=1), n.ahead=1), tail(pred_vec_ref$pred,1), check=FALSE)

############################
########### NLAR
############################

mod_ar <- linear(lynx[1:100], m=1)

predict_rolling.nlar(mod_ar, newdata=lynx[101:114])$pred[1,1]
predict_rolling.nlar(mod_ar, newdata=lynx[101:114], n.ahead=2)$pred[2,1]
predict_rolling.nlar(mod_ar, newdata=lynx[101:114], n.ahead=5)$pred[5,1]
predict_rolling.nlar(mod_ar, newdata=lynx[101:114], n.ahead=1:5)$pred[c(1,14+2,4*14+5),]

predict(mod_ar, n.ahead=5)


mod_set <- setar(lynx[1:100], m=1)
predict_rolling.nlar(mod_set, newdata=lynx[101:114])$pred[1,1]
predict_rolling.nlar(mod_set, newdata=lynx[101:114], n.ahead=5)$pred [5,1]
predict(mod_set, n.ahead=5)

predict_rolling(mod_set, newdata=lynx[101:114], n.ahead=1:3)$pred[c(1,14,15, 29),]

}