#Author: Antonio, Fabio Di Narzo. Last Modified 30 March 2011
autotriples <- function(x, lags=1:2, h, type=c("levels","persp","image", "lines", "points"), GUI=interactive()) {
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
  if(GUI) {
    require(tcltk) || stop("tcltk package is required for displaying the GUI")
    replot <- function(...) {
      type <- as.character(tclObj(ttype))
      h.new <- as.numeric(tclObj(th))
      lag1.new <- as.numeric(tclObj(tlag1))
      lag2.new <- as.numeric(tclObj(tlag2))
      if(!is.na(h.new)) h <<- h.new
      if(!is.na(h.new))
        h <<- h.new
      else
        tclvalue(th) <- as.character(h)
      if((!is.na(lag1.new))&(!is.na(lag2.new)))
        lags <<- c(lag1.new, lag2.new)
      else {
        tclvalue(tlag1) <- as.character(lags[1])
        tclvalue(tlag2) <- as.character(lags[2])
      }
      X <<- embedd(x, lags=c(-lags, 0))
      xlab <<- paste("lag",lags[1])
      ylab <<- paste("lag",lags[2])
      mod <- sm.regression(X[,1:2], X[,3], h=rep(exp(h),2), display="none")
      panel[[type]](mod$estimate)
    }
    types <- names(panel)
    ttype <- tclVar(type)
    tlag1 <- tclVar(lags[1])
    tlag2 <- tclVar(lags[2])
    th <- tclVar(h)
    frLeft <- Frame(opts=list(anchor="w"), conf=list(borderwidth=4, relief="raised"))
    frRight <- Frame(opts=list(anchor="w"), conf=list(borderwidth=4, relief="raised"))
    frRoot <- Frame(opts=list(side="left"))
    add(frRoot, frLeft, frRight)

    add(frLeft, Widget(opts=list(type="label", text="plot type:")))
    for(i in 1:length(types))
      add( frLeft, Widget(opts=list(type="radiobutton", command=replot, 
                            text=types[i], value=types[i], variable=ttype)) )

    lagEntry1 <- Widget(opts=list(type="entry", textvariable=tlag1, width=2))
    lagEntry2 <- Widget(opts=list(type="entry", textvariable=tlag2, width=2))
    add(frRight, Widget(opts=list(type="label", text="first lag:")), lagEntry1,
        Widget(opts=list(type="label", text="second lag:")), lagEntry2)

    h.start <- log(h) - 2
    h.end <- log(h) + 2
    hScale <- Widget(opts=list(type="scale", command=replot, from=h.start, to=h.end, 
                       variable=th, resolution=(h.start-h.end)/100, orient="horiz"))
    add(frRight, Widget(opts=list(type="label", text="kernel window:")), hScale)

    buildDialog(title="plot options", frRoot)
    tkbind(lagEntry1$tkvar, "<Return>",replot)
    tkbind(lagEntry2$tkvar, "<Return>",replot)
    return(invisible(NULL))
  }
  invisible(NULL)
}
