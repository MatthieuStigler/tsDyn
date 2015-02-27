
#' @S3method predict VAR
predict.VAR <- function(object, newdata, n.ahead=5, ...){
  lag <- object$lag
  k <- object$k
  include <- object$include
  if(attr(object, "varsLevel")=="ADF") stop("Does not work with VAR in diff specification")
  
  ## get coefs
  B <- coef(object)
  if(attr(object, "varsLevel")=="diff") {
    B <- VARrep.VAR(object)
    lag <- lag+1
  }
  
  ## setup starting values (data in y), innovations (0)
  original.data <- object$model[,1:k, drop=FALSE]
  starting <-   myTail(original.data,lag)
  innov <- matrix(0, nrow=n.ahead, ncol=k)
  if(!missing(newdata)) {
    if(!inherits(newdata, c("data.frame", "matrix","zoo", "ts"))) stop("Arg 'newdata' should be of class data.frame, matrix, zoo or ts")
    if(nrow(newdata)!=lag) stop("Please provide newdata with nrow=lag")
    starting <-  newdata 
  }
  
  ## use VAR sim
  res <- VAR.gen(B=B, lag=lag, n=n.ahead, 
                 starting=starting, innov=innov,include=include)
  
  ## results
  colnames(res) <- colnames(original.data )
  res <- tail(res, n.ahead)
  
  rownames(res) <- (nrow(original.data)+1):(nrow(original.data)+n.ahead)
  
  return(res)
}



if(FALSE){
  ### Predict
  environment(predict.VAR) <- environment(VECM)
  #   data(Canada)
  n <- nrow(Canada)
  
  Var_1 <- lineVar(Canada, lag=1)
  all.equal(predict(Var_1),predict(Var_1, newdata=Canada[n,,drop=FALSE]))
  all.equal(tail(fitted(Var_1),1),
            predict2(Var_1, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)
  predict2(Var_1, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE])
  Var_2 <- lineVar(Canada, lag=2)
  all.equal(predict(Var_2),predict(Var_2, newdata=Canada[c(n-1,n),,drop=FALSE]))
  all.equal(tail(fitted(Var_2),1),predict(Var_2, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)
  
}