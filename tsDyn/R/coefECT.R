#' Extract cointegration parameters A, B and PI
#' 
#' Extract parameters in VECM: adjustment coefficients \code{A}, 
#' cointegrating coefficients \code{B} , or the composite matrix \code{PI}
#' @param object An object of class \code{\link{VECM}}, \code{\link[urca]{ca.jo}}
#' @param \ldots Further arguments passed to methods
#' @author Matthieu Stigler
#' @details The functions extract the parameters from a VECM with \eqn{K} variables 
#' and rank \eqn{r}:
#'   \describe{
#'    \item{A}{Adjustment coefficients, of dim \eqn{K \times r}}
#'    \item{B}{Cointegrating coefficients, of dim \eqn{K \times r}}
#'    \item{Pi}{Matrix \eqn{\Pi=A B^{'}}, of dim \eqn{K \times K}}
#'    }
#'    Coefficients are extracted from a VECM in package \code{tsDyn}, or from a VECM 
#'    obtained in package \code{urca} from \code{\link[urca]{ca.jo}} or \code{\link[urca]{cajorls}}. 
#'    
#'    Note that by default, the A and B coefficients returned are normalized (see below). This is 
#'    the case for results obtained from \code{\link{VECM}}/\code{\link{lineVar}} and 
#'    \code{\link[urca]{cajorls}}, while for \code{\link[urca]{ca.jo}}, the user has the choice
#'    (but normalize=TRUE by default), in which case the rank \code{r} is also to be specified.
#'    The normalization is the Phillips triangular representation,  as suggested by Johansen (1995, p. 72), 
#'    standardising  the first \eqn{r\times r} coefficients to \eqn{I_r}:
#'    \describe{
#'      \item{B}{\eqn{B_{norm}=B (c^{'}B)^{-1}} with \eqn{c=(I_r,0_{K-r,r})^{'}}}
#'      \item{A}{\eqn{A_{norm}=AB^{'}c}}
#'    }
#'    Finally, note that the function also apply to objects obtained from tests of class 
#'    \code{ca.jo.test} (from \code{\link[urca]{blrtest}} etc...). Care should be taken 
#'    however, since the normalization might override the restrictions imposed. 
#' @return A matrix containing the coefficients 
#' @references Johansen, Soren, (1995), Likelihood-Based Inference in Cointegrated Vector Autoregressive Models, Oxford University Press
#' @export
#' @examples
#' data(barry)
#' vecm <- VECM(barry,  lag=1, estim="ML")
#' vecm_r2 <- VECM(barry,  lag=1, estim="ML", r=2)
#' 
#' ## extract coefficients:
#' coefA(vecm)
#' coefB(vecm)
#' coefPI(vecm)
#' coefB(vecm_r2)
#' coefPI(vecm_r2)
#' 
#' ## Beta-Restricted VECM:
#' beta_vecm2 <- coefB(vecm_r2) 
#' beta_vecm2[3,2] <- 0.02
#' vecm_r2_rest <- VECM(barry,  lag=1, estim="ML", r=2, beta=beta_vecm2)
#' round(coefB(vecm_r2_rest),5)
#' 
#' ## Package vars/urca
#' if(require(urca)){
#'  vecm_ur <- ca.jo(barry, K=2)
#'  coefB(vecm_ur)
#'  coefB(vecm_ur,r=2)
#'  coefB(cajorls(vecm_ur, r=2))
#'  all.equal(coefB(vecm), coefB(vecm_ur), check.attributes=FALSE)
#'  all.equal(coefB(vecm_r2), coefB(vecm_ur, r=2), check.attributes=FALSE)
#' }


coefB <- function(object, ...) UseMethod("coefB")

#' @rdname coefB
#' @export
coefB.VECM <- function(object,...){
  object$model.specific$beta
}

#' @export
coefB.list <- function(object,...){
  if(all(names(object)==c("rlm","beta"))){
    return(object$beta)
  } else {
    stop("No method for this object")
  }
}

#' @rdname coefB
#' @export
# Normalize by b (Z')-1, with  (Z')-1 = (C'B)-1, c: (diag|0) 
#' @param r The cointegrating rank
#' @param normalize Whether to normalize the A/B coefficients. See details
coefB.ca.jo <- function(object,r=1, normalize=TRUE, ...){
  beta <- object@V[,1:r,drop=FALSE]
  if(r>1&& normalize){
    C1 <- diag(r)
    C2 <- matrix(0, nrow = nrow(beta) - r, ncol = r)
    C <- rbind(C1, C2)
    beta <- beta %*% solve(t(C) %*% beta)
  }
  beta
}

#' @export
coefB.cajo.test <- function(object,r=1, normalize=TRUE, ...) 
  coefB.ca.jo(object=object,r=r, normalize=normalize, ...)

#' @rdname coefB
#' @export
coefA <- function(object, ...) UseMethod("coefA")

#' @rdname coefB
#' @export
coefA.VECM <- function(object,...){
  r <- object$model.specific$r
  coef(object)[,1:r, drop=FALSE]
}

#' @export
coefA.list <- function(object,...){
  if(all(names(object)==c("rlm","beta"))){
    r <- ncol(object$beta)
    return(t(coef(object$rlm)[1:r,,drop=FALSE]))
  } else {
    stop("No method for this object")
  }
}

#' @rdname coefB
#' @export
# BETA: Normalize by b (Z')-1, with  (Z')-1 = (C'B)-1, c: (diag|0)
# ALPHA: Normalize by Z, i.e. B'C, 
coefA.ca.jo <- function(object,r=1, normalize=TRUE, ...){
  alpha <- object@W[,1:r,drop=FALSE]
  if(r>1&& normalize){
    beta <- object@V[,1:r,drop=FALSE]
    C1 <- diag(r)
    C2 <- matrix(0, nrow = nrow(alpha) - r, ncol = r)
    C <- rbind(C1, C2)
    alpha <- alpha %*% t(beta)  %*% C
  }
  alpha
}

#' @export
coefA.cajo.test <- function(object,r=1, normalize=TRUE, ...) 
  coefA.ca.jo(object=object,r=r, normalize=normalize, ...)

#' @rdname coefB
#' @export
coefPI <- function(object, ...) UseMethod("coefPI")

#' @export
coefPI.default <- function(object, ...){
  coefA(object)%*%t(coefB(object))
}

#' @export
coefPI.ca.jo <- function(object,r=1, normalize=TRUE, ...){
  alpha <- object@W[,1:r,drop=FALSE]
  beta <- object@V[,1:r,drop=FALSE]
  alpha%*%t(beta)
}    

#' @export
coefPI.cajo.test <- function(object,r=1, normalize=TRUE, ...) 
  coefPI.ca.jo(object=object,r=r, normalize=normalize, ...)

if(FALSE){
  #library(vars)
  #library(tsDyn)
  #data(barry)
  H1 <- ca.jo(barry, type='trace', K=2)
  H1_ts <- VECM(barry,  lag=1, estim="ML")
  H1_ts_r2 <- VECM(barry,  lag=1, estim="ML", r=2)
  beta_r <- coefB(H1_ts_r2) 
  beta_r[3,2] <- 0.01
  H1_ts_r2_res <- VECM(barry,  lag=1, estim="ML", r=2, beta=beta_r)
  all(coefB(H1_ts_r2_res)==beta_r)
  coefA(H1_ts_r2_res)
  
  ## test coefB
  
  # r=1
  all.equal(coefB(H1_ts), coefB(cajorls(H1)), check.attributes=FALSE)
  all.equal(coefB(H1_ts), coefB(H1), check.attributes=FALSE)

  # r=2
  all.equal(coefB(H1_ts_r2), coefB(cajorls(H1,r=2)), check.attributes=FALSE)
  all.equal(coefB(H1_ts_r2), coefB(H1,r=2), check.attributes=FALSE)
  
  ## coefA
  all.equal(coefA(H1_ts),  coefA(H1, r=1) , check.attributes=FALSE)
  all.equal(coefA(H1_ts),  coefA(cajorls(H1,r=1)) , check.attributes=FALSE)
  all.equal(coefA(H1_ts_r2),  coefA(H1, r=2) , check.attributes=FALSE)
  all.equal(coefA(H1_ts_r2),  coefA(cajorls(H1,r=2)) , check.attributes=FALSE)
 
  # CoefPI
  all.equal(coefPI(H1_ts), coefPI(cajorls(H1, r=1)), check.attributes=FALSE)
  all.equal(coefPI(H1_ts), coefPI(H1, r=1), check.attributes=FALSE)
  all.equal(coefPI(H1_ts_r2), coefPI(cajorls(H1, r=2)), check.attributes=FALSE)
  all.equal(coefPI(H1_ts_r2), coefPI(H1, r=2), check.attributes=FALSE)
}

