extend <- function(...)
	UseMethod("extend")

extend.list <- function(this, subclass, ..., listV) {
	class(this) <- c(subclass, class(this))
	if(missing(listV))
		return(structure(c(this, list(...)), class=class(this)))
	else
		return(structure(c(this, listV), class=class(this)))
}

MyEnv <- function(...) {
	this <- new.env()
	vars <- list(...)
	nms <- names(vars)
	if(length(nms)) for(i in 1:length(nms))
		assign(nms[i], vars[[i]], env=this)
	structure(this, class=c("MyEnv",class(this)))
}

extend.MyEnv <- function(this, subclass, ...) {
	class(this) <- c(subclass, class(this))
	newvars <- list(...)
	if(length(names(newvars))) for(nm in names(newvars))
		assign(nm, newvars[[nm]], env=this)
	return(this)
}

availableModels <- function()
	fitters

latex <- function(obj, ...)
	UseMethod("latex")

formatSignedNum <- function(x, ...) {
  signChar <- c("-","+")[(x>=0)+1]
  nm <- names(x)
  res <- paste(signChar,format(abs(x), ...) )
  names(res) <- nm
  return(res)
}

build <- function(...)
	UseMethod("build")

add <- function(...)
	UseMethod("add")

sigmoid <- function(x) 1/(1 + exp(-x))

dsigmoid <- function(x) sigmoid(x) * (1 - sigmoid(x))

d2sigmoid <- function(x) dsigmoid(x) * (1 - 2 * sigmoid(x))

repmat <- function(a, b, c) kronecker(matrix(1,b,c), a)
