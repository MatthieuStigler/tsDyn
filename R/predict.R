

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
      res[i] <- tsDyn:::oneStep(object, newdata = t(as.matrix(res[i - xrange])),
			itime=(i - n.used), ...)
      res[i] <- res[i] + innov[i-n.used]
    }
    return(res)
  }

  res <- replicate(nboot, pred_fun(res=res))
  res_means <- rowMeans(res, na.rm=TRUE)
  pred <- res_means[n.used + 1:n.ahead]

## Compute SE if MC/boot
  if(type!="naive"){
    SE <- t(apply(res[n.used + 1:n.ahead, ,drop=FALSE], 1 ,quantile, prob=c(1-ci, ci), na.rm=TRUE))
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

set <- setar(lynx, m=2)
tsDyn:::predict.nlar(set, n.ahead=1)

environment(predict.nlar) <- environment(star)
environment(oneStep) <- environment(star)

## n-ahead: 1
tsDyn:::predict.nlar(set, n.ahead=1)
tsDyn:::predict.nlar(set, n.ahead=1, type="naive")
tsDyn:::predict.nlar(set, n.ahead=1, type="MC")
tsDyn:::predict.nlar(set, n.ahead=1, type="bootstrap")

## n-ahead: 2
tsDyn:::predict.nlar(set, n.ahead=2)
tsDyn:::predict.nlar(set, n.ahead=2, type="naive")
tsDyn:::predict.nlar(set, n.ahead=2, type="MC", nboot=100)
tsDyn:::predict.nlar(set, n.ahead=2, type="bootstrap", nboot=100)

## n-ahead: 5
tsDyn:::predict.nlar(set, n.ahead=5)
tsDyn:::predict.nlar(set, n.ahead=5, type="naive")
tsDyn:::predict.nlar(set, n.ahead=5, type="MC", nboot=5)

###### LSTAR
mod.lst <- lstar(lynx, m=2, th=1300)

predict(mod.lst, n.ahead=5)
predict(mod.lst, n.ahead=5, type="MC", nboot=100)
predict(mod.lst, n.ahead=5, type="bootstrap", nboot=100)

###### AAR
mod.aar <- aar(lynx, m=2)
predict(object=mod.aar , n.ahead=5)
}
