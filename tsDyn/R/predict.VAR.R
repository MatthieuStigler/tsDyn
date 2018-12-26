#' Predict method for objects of class \sQuote{\code{VAR}}, \sQuote{\code{VECM}} or \sQuote{\code{TVAR}}
#' 
#' Forecating the \strong{level} of a series estimated by  \sQuote{\code{VAR}} / \sQuote{\code{VECM}} or \sQuote{\code{TVAR}}
#' 
#' @aliases  predict.VAR predict.VECM
#' @param object An object of class  \sQuote{\code{VAR}}, \sQuote{\code{VECM}} or \sQuote{\code{TVAR}}
#' @param newdata Optional. A new data frame to predict from. 
#' This should contain lags of the level of the original series. See Details. 
#' @param n.ahead An integer specifying the number of forecast steps.
#' @param exoPred vector/matrix of predictions for the exogeneous variable(s) (with \sQuote{\code{n.ahead}} rows).
#' Only for \sQuote{\code{VAR}}/\sQuote{\code{VECM}}, not for \sQuote{\code{TVAR}}.
#' @param newdataTrendStart If \sQuote{\code{newdata}} is provided by the user, 
#' and the estimated model includes a trend, 
#' this argument specifies where the trend should start
#' @param \dots Arguments passed to the unexported \sQuote{\code{VAR.gen}} or \sQuote{\code{TVAR.gen}} function
#' 
#' @details
#' The forecasts are obtained recursively, and are for the levels of the series.  
#' 
#' When providing newdata, newdata has to be ordered chronologically, 
#' so that the first row/element is the earliest value. 
#' 
#' For VECM, the forecasts are obtained by transforming the VECM to a VAR (using function \code{\link{VARrep}}). 
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
#' data(barry)
#' barry_in <- head(barry, -5)
#' barry_out <- tail(barry, 5)
#' 
#' mod_vecm <- VECM(barry_in, lag=2)
#' mod_var <- lineVar(barry_in, lag=3)
#' mod_tvar <- TVAR(barry_in, lag=3, nthresh=1, thDelay=1)
#' 
#' pred_vecm <- predict(mod_vecm)
#' pred_var  <- predict(mod_var) 
#' pred_tvar <- predict(mod_tvar)
#' 
#'  
#' ## compare forecasts on a plot
#' n <- 30
#' plot(1:n, tail(barry[,1], n), type="l", xlim=c(0,n))
#' lines((n-5+1):n, pred_var[,1], lty=2, col=2)
#' lines((n-5+1):n, pred_vecm[,1], lty=2, col=3)
#' lines((n-5+1):n, pred_tvar[,1], lty=2, col=4) 
#' legend("bottomright", lty=c(1,2,2,2), col=1:4, legend=c("true", "var", "vecm", "tvar"))
#' 
#' ## example for newdata:
#' all.equal(predict(mod_vecm), predict(mod_vecm, newdata=barry[c(317, 318, 319),]))



###################
#' @export
predict.VAR <- function(object, newdata, n.ahead=5, 
                        newdataTrendStart, exoPred=NULL, ...){
  
  #   model=c("VAR", "VECM"), model <- match.arg(model)
  model <- ifelse(inherits(object, "VECM"), "VECM", "VAR")
  lag <- object$lag
  k <- object$k
  include <- object$include
  
  ## Check exogen
  hasExo <- object$exogen
  if(hasExo){
    if(is.null(exoPred)) stop("Please provide exogeneous values. ")
    if(!is.matrix(exoPred)) exoPred <- as.matrix(exoPred)
    dim_exo <- c(n.ahead, object$num_exogen)
    if(!all(dim(exoPred)==dim_exo)) 
      stop("exoPred should be a matrix with ", n.ahead, " rows and ", object$num_exogen, " columns")
  }
  
  
  ## get coefs
  if(model=="VAR"){
    B <- coef(object)
    if(attr(object, "varsLevel")=="ADF") stop("Does not work with VAR in diff specification")
    if(attr(object, "varsLevel")=="diff") {
      B <- VARrep.VAR(object)
      lag <- lag+1
    }  
  } else {
    B <- VARrep(object)
    lag <- lag+1
    
    ## check deterministc specification
    LRinclude <- object$model.specific$LRinclude
    if(LRinclude!="none"){
      if(LRinclude=="const"){
        include <- "const"
      } else if(LRinclude=="trend"){
        include <- if(include=="const") "both" else "trend"
      } else if(LRinclude=="both"){
        include <- "both"
      }
    }
    
  }
  
  
  ## setup starting values (data in y), innovations (0)
  original.data <- object$model[,1:k, drop=FALSE]
  starting <-   if(lag>0) myTail(original.data,lag) else NULL
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
  res <- VAR.gen(B=B, lag=lag, n=n.ahead, trendStart=trendStart,
                 starting=starting, innov=innov,include=include,
                 exogen=exoPred, ...)
  
  ## results
  colnames(res) <- colnames(original.data )
  
  end_rows <- nrow(original.data) + n.ahead
  if(hasArg("returnStarting") && isTRUE(list(...)["returnStarting"])) {
    start_rows <- nrow(original.data)+1 - lag
  } else {
    start_rows <- nrow(original.data)+1
  }
  rownames(res) <- start_rows : end_rows
  return(res)
}

# predict.VECM <- function(object, newdata, n.ahead=5, 
#                             newdataTrendStart, exoPred=NULL, ...){
#   
#   predict.VARVECM(object=object, newdata=newdata, n.ahead=n.ahead, 
#                   newdataTrendStart=newdataTrendStart, exoPred=exoPred,
#                   model="VECM")
#     
# }
# 
# predict.VAR <- function(object, newdata, n.ahead=5, 
#                         newdataTrendStart, exoPred=NULL, ...){
#   
#   predict.VARVECM(object=object, newdata=newdata, n.ahead=n.ahead, 
#                   newdataTrendStart=newdataTrendStart, exoPred=exoPred,
#                   model="VAR")
#   
# }



################################################
############ OLD
################################################
# unused so far, make sure;
predict2.VAR <- function(object, newdata, n.ahead=5, newdataTrendStart, exoPred=NULL, ...){
  lag <- object$lag
  k <- object$k
  include <- object$include
  hasExo <- object$exogen
  if(hasExo&&is.null(exoPred)) stop("Please provide exogeneous values. ")
  
  if(attr(object, "varsLevel")=="ADF") stop("Does not work with VAR in diff specification")
  
  ## get coefs
  B <- coef(object)
  if(attr(object, "varsLevel")=="diff") {
    B <- VARrep.VAR(object)
    lag <- lag+1
  }
  
  ## setup starting values (data in y), innovations (0)
  original.data <- object$model[,1:k, drop=FALSE]
  starting <-   if(lag>0) myTail(original.data,lag) else NULL
  innov <- matrix(0, nrow=n.ahead, ncol=k)  
  
  
  if(!missing(newdata)) {
    if(!inherits(newdata, c("data.frame", "matrix","zoo", "ts"))) stop("Arg 'newdata' should be of class data.frame, matrix, zoo or ts")
    if(nrow(newdata)!=lag) stop("Please provide newdata with nrow=lag")
    starting <-  newdata 
  }
  
  ## trend
  if(missing(newdataTrendStart)){
     if(include%in%c("trend", "both")){
       trendStart <- tail(object$model[,"Trend"],1)+1
     }  else {
       trendStart <- 0
     }
  } else {
    trendStart <- newdataTrendStart
  }
  
  
  ## use VAR sim
  res <- VAR.gen(B=B, lag=lag, n=n.ahead, trendStart=trendStart,
                 starting=starting, innov=innov,include=include,
                 exogen=exoPred, ...)
  
  ## results
  colnames(res) <- colnames(original.data )
  #   res <- tail(res, n.ahead)
  
  rownames(res) <- (nrow(original.data)+1):(nrow(original.data)+n.ahead)
  return(res)
}



if(FALSE){
  ### Predict
  environment(predict.VAR) <- environment(VECM)
  #   data(Canada)
  # data(barry)
  n <- nrow(Canada)
  
  Var_1 <- lineVar(Canada, lag=1)
  all.equal(predict(Var_1),predict(Var_1, newdata=Canada[n,,drop=FALSE]))
  all.equal(tail(fitted(Var_1),1),
            predict2(Var_1, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)
  predict2(Var_1, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE])
  
  Var_2 <- lineVar(Canada, lag=2)
  all.equal(predict(Var_2),predict(Var_2, newdata=Canada[c(n-1,n),,drop=FALSE]))
  all.equal(tail(fitted(Var_2),1),predict(Var_2, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)
  
  ## with trend
  barry_df <- as.data.frame(barry)
  n_bar <- nrow(barry)
  
  var_l1_co <-lineVar(barry, lag=1, include="const")
  var_l2_co <-lineVar(barry, lag=2, include="const")
  
  var_l1_tr <-lineVar(barry, lag=1, include="trend")
  var_l1_bo <-lineVar(barry, lag=1, include="both")
  
  tail(fitted(var_l1_co),2)
            predict.VAR(var_l1_co, newdata=barry[n_bar-1,,drop=FALSE], n.ahead=1)
  check.pred <- function(x,...){
    true <- tail(fitted(x),1)
    newD <- barry[nrow(barry)-(x$lag:1),,drop=FALSE]
    check <- predict.VAR(x, newdata=newD, n.ahead=1, ...)
    if(isTRUE(all.equal(true, check, check.attributes=FALSE))){
      res <- TRUE
    }else {
      res<-rbind(true, check)
      rownames(res) <-paste(c("true", "pred"), each=nrow(true))
    }
    return(res)
  }
  check.pred(var_l1_co)
  check.pred(x=var_l1_bo,trendStart=322)
  check.pred(x=var_l2_co)
  
  all.equal(,
            predict.VAR(var_l1_co, newdata=barry[n_bar-1,,drop=FALSE], n.ahead=1), 
            check.attributes=FALSE)
  
  
  all.equal(predict(var_l1_tr),predict(var_l1_tr, newdata=Canada[c(n_bar-1,n_bar),,drop=FALSE]))
  all.equal(predict(var_l1_bo),predict(var_l1_bo, newdata=Canada[c(n_bar-1,n_bar),,drop=FALSE]))
  
  ## exo:
  var_l1_coAsExo <-lineVar(barry, lag=1, include="none", exogen=rep(1, nrow(barry)))
  var_l1 <-lineVar(barry, lag=1, include="const")
  environment(predict.VAR) <- environment(VECM)
  environment(VAR.gen) <- environment(TVECM)
  all.equal(coef(var_l1_coAsExo), coef(var_l1)[, c(2:4,1)], check.attributes=FALSE)
  all.equal(predict.VAR(var_l1_coAsExo, exoPred=rep(1,5), n.ahead=5),   predict.VAR(var_l1))
}
















predict3.VECM <- function(object, newdata, n.ahead=5, newdataTrendStart, exoPred=NULL, ...){
  predict.VAR(object=object, newdata=newdata, n.ahead=n.ahead, 
                   newdataTrendStart=newdataTrendStart, exoPred=exoPred, model="VECM",...)
    
}

predict2.VECM <- function(object, newdata, n.ahead=5, newdataTrendStart, exoPred=NULL, ...){
  
  lag <- object$lag
  k <- object$k
  include <- object$include
  LRinclude <- object$model.specific$LRinclude
  hasExo <- object$exogen
  if(hasExo&&is.null(exoPred)) stop("Please provide exogeneous values. ")
  
  if(attr(object, "varsLevel")=="ADF") stop("Does not work with VAR in diff specification")
  
  ## get VAR representation
  B <- VARrep(object)
  lag <- lag+1
  
  ## check deterministc specification
  if(LRinclude!="none"){
    if(LRinclude=="const"){
      include <- "const"
    } else if(LRinclude=="trend"){
      include <- if(include=="const") "both" else "trend"
    } else if(LRinclude=="both"){
      include <- "both"
    }
  }
  
  ## setup starting values (data in y), innovations (0)
  original.data <- object$model[,1:k, drop=FALSE]
  starting <-   if(lag>0) myTail(original.data,lag) else NULL 
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
  res <- VAR.gen(B=B, lag=lag, n=n.ahead, trendStart=trendStart,
                 starting=starting, innov=innov,include=include,
                 exogen=exoPred,...)
  
  ## results
  colnames(res) <- colnames(original.data )
  
  rownames(res) <- (nrow(original.data)+1):(nrow(original.data)+n.ahead)
  return(res)
  
}

