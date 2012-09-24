

#### Predict nlar
predict.nlar <- function(object, newdata, n.ahead=1, type=c("naive", "MC", "bootstrap"), nboot=100, ci=0.95, ...)
{

  type <- match.arg(type)

  if(missing(newdata)) 
    newdata <- object$str$x

  res <- newdata
  n.used <- length(res)
  m <- object$str$m
  d <- object$str$d
  sd <- sqrt( mse(object) )
  steps <- object$str$steps
  resid <- as.numeric(residuals(object))
  resid <- resid[-is.na(resid)]

  tsp(res) <- NULL
  class(res) <- NULL

  res <- c(res, rep(0, n.ahead))
  xrange <- (m-1)*d + steps - ((m-1):0)*d

  if(type=="naive") nboot <-1

  ## Loop for new predictions using specific oneStep()
  pred_fun <- function(res){ 
    innov <- if(type=="naive") rep(0, n.ahead) else if (type=="MC") rnorm(n.ahead, mean=0, sd=sd) else sample(resid, size=n.ahead, replace=TRUE)
    for(i in n.used + (1:n.ahead)) {
      res[i] <- tsDyn:::oneStep.setar(object, newdata = t(as.matrix(res[i - xrange])),
			itime=(i - n.used), ...)
      res[i] <- res[i] + innov[i-n.used]
    }
    return(res)
  }

  res <- replicate(nboot, pred_fun(res=res))
  res_means <- rowMeans(res)
  pred <- res_means[n.used + 1:n.ahead]

## Compute SE if MC/boot
  if(type!="naive"){
    SE <- t(apply(res[n.used + 1:n.ahead, ,drop=FALSE], 1 ,quantile, prob=c(1-ci, ci)))
    SE <- ts(SE, start = tsp(newdata)[2] + deltat(newdata),
             frequency=frequency(newdata))
  }

## Return result
  pred <- ts(pred, start = tsp(newdata)[2] + deltat(newdata),
             frequency=frequency(newdata))
  if(type=="naive"){
    result <- pred
  } else {
    result <- list(pred=pred, se=SE)
  }

  return(result)
  
}

if(FALSE){
library(tsDyn)

set <- setar(lynx, m=1)
tsDyn:::predict.nlar(set, n.ahead=1)

environment(predict2.nlar) <- environment(star)
environment(oneStep) <- environment(star)

## n-ahead: 1
tsDyn:::predict.nlar(set, n.ahead=1)
predict2.nlar(set, n.ahead=1, type="naive")
predict2.nlar(set, n.ahead=1, type="MC")
predict2.nlar(set, n.ahead=1, type="bootstrap")

## n-ahead: 2
tsDyn:::predict.nlar(set, n.ahead=2)
predict2.nlar(set, n.ahead=2, type="naive")
predict2.nlar(set, n.ahead=2, type="MC", nboot=100)
predict2.nlar(set, n.ahead=2, type="bootstrap", nboot=100)

## n-ahead: 5
tsDyn:::predict.nlar(set, n.ahead=5)
predict2.nlar(set, n.ahead=5, type="naive")
predict2.nlar(set, n.ahead=5, type="MC", nboot=5)

}