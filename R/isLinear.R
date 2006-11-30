isLinear <- function(object, ...)
  UseMethod("isLinear")

isLinear.default <- function(object, ...)
  stop("no linearity tests available for this model")

