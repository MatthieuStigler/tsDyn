vec2var.tsDyn <- function(vecm) {

  lag <- vecm$lag
  k <- vecm$k
  r <- vecm$model.specific$r
  co <- vecm$coefficients

##obtain matrices
  ## Pi matrix
  Pi <-  co[, grep("ECT", colnames(co))]%*%t(vecm$model.specific$beta)
  ## A_i matrix
  Amat <- matrix(NA, nrow=k, ncol=lag*(k)+k)
  ## A_lag+1 matrix
  Amat[,(1:k)+k*lag] <- -co[,grep(paste("-", lag, sep=""), colnames(co))]
  ## A_lag+1 matrix
  if(lag>1) for(i in 1:(lag-1)) Amat[,(1:k)+k*i] <- -(co[,grep(paste("-", i, sep=""), colnames(co))] -co[,grep(paste("-", i+1, sep=""), colnames(co))])

  cumulMat <- matrix(0, k,k)
  for(i in 1:lag) cumulMat <- cumulMat + Amat[,(1:k)+k*i]
  Amat[, 1:k] <- Pi + (diag(k)- cumulMat )

## Names
  rownames(Amat) <- rownames(co)
  varNames <- colnames(vecm$model)[i:k]
  colnames(Amat) <- paste(rep(varNames, k), rep(1:(lag+1), each=k), sep=" -")

## Add deterministic terms
  if(vecm$include!="none"){
    incName <- switch(vecm$include, "const"="Intercept", trend="Trend", both="Intercept|Trend")
    Amat <- cbind(co[,grep(incName , colnames(co)),drop=FALSE],Amat)
  }

## res
  Amat
}


########### EXAMPLE
if(FALSE){

data(denmark)
dat_examp <- denmark[,2:3]

vec2var.tsDyn(VECM(dat_examp, lag=1, include="const", estim="ML"))
  vec1 <- vec2var(ca.jo(dat_examp,  K=2, spec="transitory"))$A
  cbind(vec1$A1, vec1$A2)

vec2var.tsDyn(VECM(dat_examp, lag=2, include="const", estim="ML"))
vec2 <- vec2var(ca.jo(dat_examp,  K=3, spec="transitory"))$A
  cbind(vec2$A1, vec2$A2, vec2$A3)




}