## Copyright (C) 2005/2006  Antonio, Fabio Di Narzo
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

delta <- function(x, m, d=1, eps) {
	if(m<2)
          stop("embedding dimension 'm' should be at least equal to 2")
	C <- d2(x, m=m+1, d=d, t=1, eps.min=eps)[1, 2:(m+2)]		#computed sample correlation integrals
	return( 1 - ( (C[m])^2 /  (C[m-1]*C[m+1]) ) )
}

delta.test <- function(x, m=2:3, d=1, eps = seq(0.5*sd(x),2*sd(x),length=4), B=49) {
	delta.b <- numeric(B)
	p.value <- matrix(NA,length(m),length(eps))
        for(j in 1:length(m)) for(i in 1:length(eps)){
		delta.c <- delta(x, m=m[j], d=d, eps=eps[i])
		for(b in 1:B)
                  delta.b[b] <- delta(sample(x), m=m[j], d=d, eps=eps[i])
		p.value[j,i] <- (1+sum(delta.b>=delta.c))/(1+B)
        }
	dimnames(p.value) <- list(m=m,eps=format(eps, digits=4))
	return(p.value)
}

delta.lin <- function(x, m, d=1) {
	V1 <- var(embedd(x, m=m+1, d=d))
	V2 <- var(embedd(x, m=m, d=d))
	tmp <- eigen(V1, symmetric=TRUE)$values[1] / eigen(V2,symmetric=TRUE)$values[1]
	return(1-tmp)
}

delta.lin.test <- function(x, m=2:3, d=1, eps = seq(0.5*sd(x),2*sd(x),length=4), B=49) {
	mu <- function(x, m, eps)
          delta(x, m=m, d=d, eps=eps) - delta.lin(x, m=m, d=d)
	mu.c <- numeric()
	n <- length(x)
	ar.model <- ar(x)	#automatic AR order selection based on AIC
	mu.b <- numeric(B)
	p.value <- matrix(NA,length(m),length(eps))
        for(j in 1:length(m)) for(i in 1:length(eps)) {
          mu.c[i] <- mu(x, m[j], eps[i])
          for(b in 1:B) {
            xb <- arima.sim(n = n, list(ar=ar.model$ar),
                            rand.gen = function(n, ...) rnorm(n, 0, sqrt(ar.model$var.pred) ) )
            mu.b[b] <- mu(xb, m[j], eps[i])
          }
          p.value[j,i] <- (1+sum(mu.b>=mu.c[i]))/(1+B)
	}
	dimnames(p.value) <- list(m=m, eps=format(eps, digits=4))
	return(p.value)
}
