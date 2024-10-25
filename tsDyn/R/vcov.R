
#' @export
vcov.TVAR <- function (object, ...) {
  x <- object
  Z<-t(as.matrix(tail.matrix(x$model[,-c(1:x$k)],x$t)))
  R <- chol2inv(chol(tcrossprod(Z)))
  
  ## 
  Sigma <- matrix((1/x$df.residual) * crossprod(x$residuals),ncol=x$k)
  VarCovB <- Sigma %x% R
  
  VarCovB
}


