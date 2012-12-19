#Author: Antonio, Fabio Di Narzo. Last Modified 30 March 2011
autopairs <- function(x, lag=1, h,
                      type=c("levels","persp","image","lines","points","regression")) {
  panel <- list(levels = function()  sm.density(X, h=rep(h,2), xlab=xlab, ylab=ylab, main="density", display="slice"),
		persp = function() sm.density(X, h=rep(h,2), xlab=xlab, ylab=ylab, main="density", display="persp"),
		image = function() sm.density(X, h=rep(h,2), xlab=xlab, ylab=ylab, main="density", display="image"),
		lines = function() plot(X, xlab=xlab, ylab=ylab, main="lines", type="l"),
		points = function() plot(X, xlab=xlab, ylab=ylab, main="scatter"),
		regression = function() sm.regression(X[,1], X[,2], h=h, xlab=xlab, ylab=ylab, main="regression", ask=FALSE))
  require(sm) || stop("sm package is required for kernel estimations")
  lags <- c(-lag, 0)
  X <- embedd(x, lags=lags)
  xlab <- paste("lag",lag)
  ylab <- paste("lag",0)
  type <- match.arg(type)
  if(missing(h)) {
    h <- hnorm(X)[1]
  }
  panel[[type]]()
}
