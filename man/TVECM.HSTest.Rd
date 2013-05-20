\name{TVECM.HStest }
\alias{TVECM.HStest }
\title{Test of linear cointegration vs threshold cointegration}
\description{
Tests the null of linear cointegration against threshold cointegration following Hansen and Seo (2002). Fixed regressor anfd residual bootstrap are available. 
}
\usage{
TVECM.HStest(data, lag=1, ngridTh=300, trim=0.05, 
 nboot=100, fixed.beta=NULL,  intercept=TRUE, 
 boot.type=c("FixedReg", "ResBoot"), 
 hpc=c("none", "foreach"))
}

\arguments{
\item{data}{Time series }
\item{lag}{Number of lags to include in each regime}
\item{ngridTh}{Number of threshold points to estimate}
\item{trim}{Trimming parameter indicating the minimal percentage of observations in each regime}
\item{nboot}{Number of bootstrap replications}
\item{fixed.beta}{Numeric. User pre-specified cointegrating value. When NULL (default), the value is estimated from the linear VECM.}
\item{intercept}{Logical. Whether an intercept has to be included in the VECM}
\item{boot.type}{Character. Type of bootstrap simulation (only if nboot>0)}
\item{hpc}{Possibility to run the bootstrap on parallel core. See details}
}

\details{
This test follows the implementation done by Hansen and Seo (2002). The cointegrating value is estimated from the linear VECM. Then, conditional on this value, the LM test is run for a range of different threshold values. The maximum of those LM test values is reported. 

Two bootstrap are available: a fixed regressor, as well as a usual residual bootstrap (using the function \code{\link{TVECM.sim}}). 

Available methods are \code{print()}, \code{summary()} and \code{plot()}. 

With the argument \code{hpc}, the burdensome bootstrap replication can be run on parallel cores, thus alleviating the time of computation. The user needs simply to choose a package (among doMC, doMPI, doSNOW or doRedis) and register the backend. See the vignette for more details. 
}

\value{A list cointaining diverse values:
\item{stat}{The sup-LM statistic.}
\item{values}{The whole LM values.}
\item{PvalBoot}{The bootstrap p-value}
\item{CriticalValBoot}{The bootstrap critical values}
\item{allBoots}{The boot sup-LM values}
\item{args}{Some user given args (nboot, boot.type)}
}
\section{Reproducibility}{Comparison with original paper is made difficult as values of the test are not shown in the paper, only their critical values, which depend on random bootstrap. 

Comparison is done with the GAUSS code available on the page of Bruce Hansen. Running tar_ci, we have the same sup-LM value when lags=1 and lags=2, a higher value with lag=3. When the test is run with pre-specified beta values, we have different results, sometimes higher but also smaller sup-LM value. 
}
\seealso{
\code{\link{zeroyld}}: data used in the original paper of Hansen and Seo. 

\code{\link{TVECM.SeoTest}}: a similar test, but with null hypothesis of no-cointegration. 

\code{\link{TVECM}} for estimating a TVECM, \code{\link{TVECM.sim}} for simulating/bootstrap a TVECM, 
}


\author{Matthieu Stigler}
\keyword{ ts }
\references{
Hansen, B. and Seo, B. (2002), Testing for two-regime threshold cointegration in vector error-correction models, Journal of Econometrics, 110, pages 293 - 318
}
\examples{
#Use original data from paper:
data(zeroyld)
dataPaper<-zeroyld
# Test: nboot, number of bootstrap replications, should be high 
\dontrun{
test1<-TVECM.HStest(dataPaper, lag=1, intercept=TRUE, nboot=1000)
}

#we use here for the example a much smaller number of bootstrap:
test1<-TVECM.HStest(dataPaper, lag=1, intercept=TRUE, nboot=10)

test1
summary(test1)
plot(test1)

#can have only specific plots:
plot(test1, which="LM values")
plot(test1, which="Density")

## Run the function in parallel:
\dontrun{
#we show here the use with package doMC
library(doMC)
registerDoMC(2) #Number of cores
test1<-TVECM.HStest(dataPaper, lag=1, intercept=TRUE, nboot=1000, hpc="foreach")
}

}


