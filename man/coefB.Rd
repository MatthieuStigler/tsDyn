\name{coefB}
\alias{coefA}
\alias{coefB}
\alias{coefPI}
\title{Extract cointegration parameter A, B and PI}
\usage{
coefB(object, ...)

coefA(object, ...)

coefPI(object, ...)
}
\arguments{
  \item{object}{An object of class \code{\link{VECM}},
  \code{\link[vars]{ca.jo}}}

  \item{\ldots}{}
}
\value{
A matrix containing the coefficients
}
\description{
Extract cointegration parameter A, B and PI
}
\examples{
data(barry)
vecm <- VECM(barry,  lag=1, estim="ML")
vecm_r2 <- VECM(barry,  lag=1, estim="ML", r=2)

## extract coefficients:
coefA(vecm)
coefB(vecm)
coefPI(vecm)
coefB(vecm_r2)
coefPI(vecm_r2)

## Beta-Restricted VECM:
beta_vecm2 <- coefB(vecm_r2)
beta_vecm2[3,2] <- 0.02
vecm_r2_rest <- VECM(barry,  lag=1, estim="ML", r=2, beta=beta_vecm2)
round(coefB(vecm_r2_rest),5)

## Package vars/urca
if(require(urca)){
 vecm_ur <- ca.jo(barry, K=2)
 coefB(vecm_ur)
 coefB(vecm_ur,r=2)
 coefB(cajorls(vecm_ur, r=2))
 all.equal(coefB(vecm), coefB(vecm_ur), check.attributes=FALSE
 all.equal(coefB(vecm_r2), coefB(vecm_ur, r=2), check.attributes=FALSE)
}
}

