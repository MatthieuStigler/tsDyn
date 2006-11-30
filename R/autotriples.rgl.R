#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2006-03-04 22:04:02 +0100 (sab, 04 mar 2006) $
autotriples.rgl <- function(x, lags=1:2, type=c("lines","points")) {
	require(rgl) || stop("rgl package is required for interactive 3d visualization")
	type <- match.arg(type)
	X <- embedd(x, lags=c(-lags,0))
	rgl.clear()
	if(type=="lines")
		rgl.linestrips(X[,1],X[,2],X[,3])
	else if (type=="points")
		rgl.points(X[,1],X[,2],X[,3])
        invisible(NULL)
}
