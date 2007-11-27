## Copyright (C) 2006  Antonio, Fabio Di Narzo
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

#FAR fitter
far <- function(str, ...){
  xx <- str$xx
  yy <- str$yy
  dat <- data.frame(cbind(xx,y=yy))
  predNames <- names(dat)
  form <- as.formula(paste("y ~", predNames))
  model <- list()
  model$internals <- loess(form, data=dat)
  model$n.used <- model$internals$n.used
  model$residuals <- model$internals$residuals
  model$fitted.values <- model$internals$fitted.values
  model$k <- model$internals$rank
  return(extend(nlar(str, ...), "far", listV=model))
}