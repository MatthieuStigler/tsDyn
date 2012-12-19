#Author: Antonio, Fabio Di Narzo. Last Modified 30 March 2011
autotriples <- function(x, lags=1:2, h, type=c("levels","persp","image", "lines", "points")) {
  require(sm) || stop("sm package is required for kernel density estimation")
  require(scatterplot3d) ||	stop("the scatterplot3d package is required for 3d visualization")
	
  panel <- list(levels = function(x) contour(x, xlab=xlab, ylab=ylab),
		persp = function(x) persp(x, xlab=xlab, ylab=ylab, zlab=zlab),
		image = function(x) image(x, xlab=xlab, ylab=ylab),
		lines = function(x) scatterplot3d(X, xlab=xlab, ylab=ylab, zlab=zlab, main="directed lines", type="l"),
		points = function(x) scatterplot3d(X, xlab=xlab, ylab=ylab, zlab=zlab, main="cloud", pch=1))
  type <- match.arg(type)
  X <- embedd(x, lags=c(-lags,0))
  if(missing(h)) 
    h <- hnorm(X[,1])
  xlab <- paste("lag",lags[1])
  ylab <- paste("lag",lags[2])
  zlab <- "lag 0"
  mod <- sm.regression(X[,1:2], X[,3], h=rep(h,2), display="none")
  panel[[type]](mod$estimate)
  invisible(NULL)
}
