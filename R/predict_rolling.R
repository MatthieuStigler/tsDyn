
predict_rolling  <- function (object, ...)  
  UseMethod("predict_rolling")

predict_rolling.default  <- function (object, ...)  NULL


predict_rolling.nlVar <- function(object, nroll=10, n.ahead=1, refit.every, newdata, ...){

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
#   if(inherits(origSerie, "ts")){
#     attributes(R) <- attributes(outSerie)
#   }

return(R)

}




##############################################################
##################### TESTS
##############################################################
if(FALSE){
library(tsDyn)
data(barry)
n_ca<- nrow(barry)

environment(predict_rolling.nlVar) <- environment(star)

#### No refit lag=1
pred_roll<-predict_rolling.nlVar(object=lineVar(barry, lag=1), nroll=10, n.ahead=1)

pred1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[n_ca-10,,drop=FALSE])
pred2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[n_ca-9,,drop=FALSE])
all.equal(rbind(pred1, pred2), head(pred_roll,2), check=FALSE)

all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[n_ca-1,,drop=FALSE]), tail(pred_roll,1), check=FALSE)

#### No refit lag=3
pred_roll<-predict_rolling.nlVar(object=lineVar(barry, lag=3), nroll=10, n.ahead=1)

pred1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-2):n_ca)-10,,drop=FALSE])
pred2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-2):n_ca)-9,,drop=FALSE])
all.equal(rbind(pred1, pred2), head(pred_roll,2), check=FALSE)

all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-2):n_ca)-1,,drop=FALSE]), tail(pred_roll,1), check=FALSE)



### Refit lag=1
pred_ref<-predict_rolling.nlVar(object=lineVar(barry, lag=1), nroll=10, n.ahead=1, refit.every=1)

pred_ref_1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1)
pred_ref_2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-9), lag=1), n.ahead=1)
all.equal(rbind(pred_ref_1, pred_ref_2), head(pred_ref,2), check=FALSE)
all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-1), lag=1), n.ahead=1), tail(pred_ref,1), check=FALSE)

### Refit lag=3
pred_ref<-predict_rolling.nlVar(object=lineVar(barry, lag=3), nroll=10, n.ahead=1, refit.every=1)

pred_ref_1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1)
pred_ref_2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-9), lag=3), n.ahead=1)
all.equal(rbind(pred_ref_1, pred_ref_2), head(pred_ref,2), check=FALSE)
all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-1), lag=3), n.ahead=1), tail(pred_ref,1), check=FALSE)


#### No refit: VAR diff,  lag=1
pred_roll_dif<-predict_rolling.nlVar(object=lineVar(barry, lag=1, I="diff"), nroll=10, n.ahead=1)

pred1_dif <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1, I="diff"), n.ahead=1, newdata=barry[(n_ca-1):n_ca-10,,drop=FALSE])
pred2_dif <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1, I="diff"), n.ahead=1, newdata=barry[(n_ca-1):n_ca-9,,drop=FALSE])
all.equal(rbind(pred1_dif, pred2_dif), head(pred_roll_dif,2), check=FALSE)

all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1, I="diff"), n.ahead=1, newdata=barry[(n_ca-1):n_ca-1,,drop=FALSE]), tail(pred_roll_dif,1), check=FALSE)


#### VECM No refit lag=1
pred_vec_roll_l1 <- predict_rolling.nlVar(object=VECM(barry, lag=1), nroll=10, n.ahead=1)

pred_vec_l1_1 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[((n_ca-1):n_ca)-10,,drop=FALSE])
pred_vec_l1_2 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[((n_ca-1):n_ca)-9,,drop=FALSE])
all.equal(rbind(pred_vec_l1_1, pred_vec_l1_2), head(pred_vec_roll_l1,2), check=FALSE)

all.equal(predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[((n_ca-1):n_ca)-1,,drop=FALSE]), tail(pred_vec_roll_l1,1), check=FALSE)



#### VECM No refit lag=3
pred_vec_roll <- predict_rolling.nlVar(object=VECM(barry, lag=3), nroll=10, n.ahead=1)

pred_vec1 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-3):n_ca)-10,,drop=FALSE])
pred_vec2 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-3):n_ca)-9,,drop=FALSE])
all.equal(rbind(pred_vec1, pred_vec2), head(pred_vec_roll,2), check=FALSE)

all.equal(predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-3):n_ca)-1,,drop=FALSE]), tail(pred_vec_roll,1), check=FALSE)



### VECM Refit lag=1
pred_vec_ref<-predict_rolling.nlVar(object=VECM(barry, lag=1), nroll=10, n.ahead=1, refit.every=1)

pred_vec_ref_1 <- predict(VECM(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1)
pred_vec_ref_2 <- predict(VECM(tsDyn:::myHead(barry,n_ca-9), lag=1), n.ahead=1)
all.equal(rbind(pred_vec_ref_1, pred_vec_ref_2), head(pred_vec_ref,2), check=FALSE)
all.equal(predict(VECM(tsDyn:::myHead(barry,n_ca-1), lag=1), n.ahead=1), tail(pred_vec_ref,1), check=FALSE)

}