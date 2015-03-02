#' Predict method for objects of class \sQuote{\code{VAR}} or \sQuote{\code{VECM}}
#' 
#' Forecating a  \sQuote{\code{VAR}} / \sQuote{\code{VECM}}
#' 
#' @aliases  predict.VAR predict.VECM
#' @param object An object of class  \sQuote{\code{VAR}} or \sQuote{\code{VECM}}
#' @param newdata Optional. A new data frame to predict from.
#' @param n.ahead An integer specifying the number of forecast steps.
#' @param exoPred vector/matrix of predictions for the exogeneous variable(s)
#' @param \dots Arguments passed to the unexported \sQuote{\code{VAR.gen}} function
#' 
#' @details
#' The forecasts are obtained recursively from the estimated model.  
#' 
#' @return A matrix of predicted values.
#' @author Matthieu Stigler
#' @seealso  \code{\link{lineVar}} and \code{\link{VECM}}.
#' 
#' A more sophisticated predict function, allowing to do sub-sample rolling
#' predictions: \code{\link{predict_rolling}}.
#' @keywords regression
#' @examples
#' 
#' data(barry)
#' mod_vecm <- VECM(barry, lag=2)
#' predict(mod_vecm)
#' 
#' # note that when providing newdata, newdata has to be ordered chronologically, 
#' # so that the first row/element is the earliest value:
#' all.equal(predict(mod_vecm), predict(mod_vecm, newdata=barry[c(322, 323, 324),]))


# unused so far, make sure;#' @S3method predict VAR
predict.VAR <- function(object, newdata, n.ahead=5, exoPred=NULL, ...){
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
  
  ## use VAR sim
  res <- VAR.gen(B=B, lag=lag, n=n.ahead, 
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