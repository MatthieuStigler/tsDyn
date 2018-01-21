#' Predict method for objects of class \sQuote{\code{TVAR}} 
#' 
#' Predicting series estimated by  \sQuote{\code{TVAR}} 
#' 
#' @aliases  predict.TVAR
#' @param object An object of class  \sQuote{\code{TVAR}}
#' @param newdata Optional. A new data frame to predict from. 
#' This should contain lags of the level of the original series. See Details. 
#' @param n.ahead An integer specifying the number of forecast steps.
#' @param newdataTrendStart If \sQuote{\code{newdata}} is provided by the user, 
#' and the estimated model includes a trend, 
#' this argument specifies where the trend should start
#' @param \dots Arguments passed to the unexported \sQuote{\code{VAR.gen}} function
#' 
#' @details
#' The forecasts are obtained recursively, and are for the levels of the series.  For VECM, the forecasts are 
#' obtained by transforming the VECM to a VAR (using function \code{\link{VARrep}}). 
#' Note that a VECM(lag=p) corresponds to a VAR(lag=p+1), so that if the user provides newdata 
#' for a VECM(lag=p), newdata should actually contain p+1 rows. 
#' 
#' @return A matrix of predicted values.
#' @author Matthieu Stigler
#' @seealso  \code{\link{lineVar}} and \code{\link{VECM}}. \code{\link{VARrep}}
#' 
#' A more sophisticated predict function, allowing to do sub-sample rolling
#' predictions: \code{\link{predict_rolling}}.
#' @keywords regression
#' @examples
#' 
#' data(zeroyld)
#' mod_tvar <- TVAR(zeroyld, lag=2, nthresh=2, thDelay=1, trim=0.1, mTh=1, plot=TRUE)
#' pred_VECM <- predict(object=mod_tvar)
#' 
#' 
#' mod_var <- lineVar(barry, lag=3)
#' pred_VAR <- predict(mod_var)
#'  
#' ## compare
#' 
#' plot(tail(barry[,1],50), type="l", xlim=c(0,60))
#' lines(51:55,pred_VAR[,1], lty=2, col=2)
#' lines(51:55,pred_VECM[,1], lty=2, col=3)
#' 
#' 
#' # note that when providing newdata, newdata has to be ordered chronologically, 
#' # so that the first row/element is the earliest value:
#' all.equal(predict(mod_vecm), predict(mod_vecm, newdata=barry[c(322, 323, 324),]))



###################
predict.TVAR <- function(object, newdata, n.ahead=5, 
                        newdataTrendStart, ...){
  
  ## extract parameters, coefs
  lag <- object$lag
  k <- object$k
  include <- object$include
  B <- object$coeffmat  
  Thresh <- getTh(object)
  nthresh <- object$model.specific$nthresh

  
  ## setup starting values (data in y), innovations (0)
  original.data <- object$model[,1:k, drop=FALSE]
  starting <-   myTail(original.data,lag)
  innov <- matrix(0, nrow=n.ahead, ncol=k)  
  
  
  if(!missing(newdata)) {
    if(!inherits(newdata, c("data.frame", "matrix","zoo", "ts"))) stop("Arg 'newdata' should be of class data.frame, matrix, zoo or ts")
    if(nrow(newdata)!=lag) stop("Please provide newdata with nrow=lag")
    starting <-  newdata 
  }
  
  ## trend
  if(missing(newdataTrendStart)){
    if(include%in%c("trend", "both")){
      trendStart <- object$t+1
    }  else {
      trendStart <- 0
    }
  } else {
    trendStart <- newdataTrendStart
  }
  
  
  ## use VAR sim
  res <- TVAR.gen(B=B, lag=lag, n=n.ahead, trendStart=trendStart,
                 starting=starting, innov=innov,include=include, 
                 Thresh=Thresh, nthresh=nthresh, returnStarting=FALSE, ...)
  
  ## results
  colnames(res) <- colnames(original.data)
  rownames(res) <- (nrow(original.data)+1):(nrow(original.data)+n.ahead)
  return(res)
}
