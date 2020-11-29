#'Extract a variable showing the regime
#'
#'This function allows to extract the indicator variable specifying the regime
#'in which the process is at time t.
#'
#'
#'@aliases regime regime.default
#'@param object object of class \code{setar} or \code{nlVar}
#'@param initVal Logical. Whether the NA initial values should be returned.
#'Default to TRUE.
#'@param timeAttr Logical. Whether the time attributes should be returned.
#'Default to TRUE.
#'@param series Optional. A numeric vector to classify according to the model. 
#'@param \dots additional arguments to \code{regime}
#'@return Time series of same attributes as input to setar.
#'@author Matthieu Stigler
#'@keywords ts
#'@export
#'@examples
#'
#'set<-setar(lynx, m=3)
#'regime(set)
#'regime(set, time=FALSE, initVal=FALSE)
#'
#'plot(regime(set))
#'

#indicator of the regime of the obs
regime <- function (object, initVal=TRUE, timeAttr=TRUE, series = NULL, ...)  
  UseMethod("regime")

#'@export
regime.default <- function(object, initVal=TRUE, timeAttr=TRUE, series = NULL, ...)
  NULL


#' @export
regime.setar <- function(object, initVal=TRUE, timeAttr=TRUE, series = NULL, ...) {
  if(!is.null(series)) {
    thDelay <- object$model.specific$thDelay
    series_lagged <- lag_manual(series, k = thDelay+1, time_attr = timeAttr)
    reg <- as.numeric(cut(series_lagged, 
                          breaks = c(-Inf,getTh(object), Inf),
                          labels = 1: (object$model.specific$nthresh+1)))
    # reg2[1:object$str$m] <- NA
    if(timeAttr) attributes(reg) <- attributes(series)
    initVal <-  TRUE
    timeAttr <-  FALSE
  } else {
    reg <- object$model.specific$regime  
  }
  
  str <- object$str
  
  if(timeAttr){
    attributes(reg) <- object$model.specific$timeAttributes
    if(!initVal) {
      reg <- window(reg, start=time(reg)[length(str$x)-length(str$yy)+1])
    }
  } else {
    if(!initVal){
      reg <- reg[-seq_len(length(str$x)-length(str$yy))]
    }
  }
  
  return(reg)
}

#' @export
regime.nlVar <- function(object, initVal=TRUE, timeAttr=TRUE, series = NULL, ...) {
  reg <- object$model.specific$regime
  
  if(timeAttr){
    attributes(reg) <- object$model.specific$timeAttributes
    if(!initVal) {
      reg <- window(reg, start=time(reg)[object$T-object$t+1])
    }
  } else {
    if(!initVal){
      reg <- reg[-c(1:(object$T-object$t))]
    }
  }
  
  return(reg)
}

#' @rdname regime
#' @param discretize logical (default TRUE) whether the series are discretized to {1,2}, 
#' or whether regime probabilities are returned. 
#' @export
regime.lstar <- function(object, initVal=TRUE, timeAttr=TRUE, series, discretize=TRUE, ...){
  
  if(!missing(series)) stop("arg 'series' not implemented for lstar")
  thVar <- object$model.specific$thVar
  str <- object$str
  
  reg <- G(z=thVar, gamma=coef(object)["gamma"], th=getTh(object))
  
  if(discretize) {
    reg <- ifelse(reg <=0.5, 1,2)
  }
  
  if(timeAttr){
    attributes(reg) <- object$model.specific$timeAttributes
    if(initVal) {
      ans <- reg
    } else {
      ans <- window(reg, start=time(reg)[length(str$x)-length(str$yy)+1])
    }
  } else {
    if(initVal){
      ans <- reg
    } else {
      ans <- reg[-c(1:(length(str$x)-length(str$yy)))]
    }
  }
  
  return(ans)
  
}

