#Author: Antonio, Fabio Di Narzo. Last Modified 30 March 2011
autopairs <- function(x, lag=1, h,
                      type=c("levels","persp","image","lines","points","regression"), GUI=interactive()) {
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
  if(GUI) {
    require(tcltk) || stop("tcltk package is required for displaying the GUI")
    replot <- function(...) {
      type <- as.character(tclObj(vType))
      h.new <- exp(as.numeric(tclObj(vH)))
      lag.new <- as.numeric(tclObj(vLag))
      if(!is.na(h.new))
        h <<- h.new
      else
        tclvalue(th) <- as.character(log(h))
      if(!is.na(lag.new))
        lag <<- lag.new
      else
        tclvalue(vLag) <- as.character(lag)
      lags <- c(-lag, 0)
      X <<- embedd(x, lags=lags)
      xlab <<- paste("lag",lag)
      ylab <<- paste("lag",0)
      panel[[type]]()
    }
    types <- c("levels","persp","image","lines","points","regression")
    vType <- tclVar(type)
    vLag <- tclVar(lag)
    vH <- tclVar(log(h))
	
    frLeft <- Frame(opts=list(anchor="w"), conf=list(borderwidth=4, relief="raised"))
    frRight <- Frame(opts=list(anchor="w"), conf=list(borderwidth=4, relief="raised"))
    frRoot <- Frame(opts=list(side="left"))
    add(frRoot, frLeft, frRight)

    add(frLeft, Widget(opts=list(type="label", text="plot type:")))
    for(i in 1:length(types))
      add( frLeft, Widget(opts=list(type="radiobutton", command=replot, 
                            text=types[i], value=types[i], variable=vType)) )

    lagEntry <- Widget(opts=list(type="entry", textvariable=vLag, width=2))
    add(frRight, Widget(opts=list(type="label", text="lag:")),lagEntry)

    h.start <- log(h) - 2
    h.end <- log(h) + 2
    F <- function (phi1, phi2, x_t, s_t) {
      noRegimes <- dim(phi1)[1]
      local <- array(0, c(noRegimes, dim(x_t)[1]))
      local[1, ] <- x_t %*% phi1[1, ]
      for (i in 2:noRegimes) local[i, ] <- (x_t %*% phi1[i, ]) * G(s_t, gamma = phi2[i - 1, 1], th = phi2[i - 1, 2])
      return(apply(local, 2, sum))
    }
    hScale <- Widget(opts=list(type="scale", command=replot, from=h.start, to=h.end, 
                       variable=vH, resolution=(h.start-h.end)/100, orient="horiz"))
    add(frRight, Widget(opts=list(type="label", text="kernel window:")),hScale)

    buildDialog(title="plot options:", frRoot)
    tkbind(lagEntry$tkvar, "<Return>", replot)
    return(invisible(NULL))
  }
}
