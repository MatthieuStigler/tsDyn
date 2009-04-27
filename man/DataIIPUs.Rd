\name{IIPUs} 
\alias{IIPUs}
\docType{data}
\title{US monthly industrial production from Hansen (1999)}
\description{
  This data, used as example in Hansen (1999), contains the US monthly industrial production.
}
\usage{data(IIPUs)}
\format{
  A monthly time series of class ts starting in january 1960 and ending in september 1997. Note that the series ends at 1997 and not 1998 as in the paper of Hansen, even if the data was taken from hi site and the graph is exactly the same. 
}
\source{
Hansen (1999) Testing for linearity, Journal of Economic Surveys, Volume 13, Number 5, December 1999 , pp. 551-576(26)
avalaible at: \url{http://www.ssc.wisc.edu/~bhansen/papers/cv.htm}
}
\examples{
data(IIPUs)
end(IIPUs) #not same date as in the paper
plot(IIPUs)#exactly same graph as in the paper
sel<-selectSETAR(IIPUs, m=16, thDelay=5, criterion="SSR", trim=0.1, plot=FALSE)
sel #R function obtains a lower SSR with another thresold
plot(sel)
setar(IIPUs, m=16, thDelay=5, trim=0.1, th=sel$th)

sel2<-selectSETAR(IIPUs, m=16, thDelay=5, criterion="SSR", trim=0.1, plot=FALSE, nthresh=2)
sel2
#all results agree
set2<-setar(IIPUs, m=16, thDelay=5, th=sel2$th, trim=0.1)
#most of the results agree, except constant in the low regime which has opposed signs!
summary(set2)

#this is obviously a error in Hansen, see:
XX<-embed(IIPUs, 17)
Y<-XX[,1]
X<-XX[,-1]
dummyDown<-ifelse(X[,6]<= -2.5, 1,0)
sum(dummyDown)
M<-cbind(1*dummyDown,X*dummyDown )
lm(Y~M-1)

setarTest(IIPUs, m=16, thDelay=5, nboot=1, check=TRUE)
#because of the discrepency. test1vs2 does not correspond, test 1vs3 conforms
setarTest(IIPUs, m=16, thDelay=5, nboot=1, check=TRUE, test="2vs3")
#test 2vs3 is also different of the version in the article (27)
}
\keyword{datasets}
\concept{Testing for linearity}

