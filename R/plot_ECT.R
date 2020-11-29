#' Plot the Error Correct Term (ECT) response
#' 
#' This plot shows how variables in a (T)VECM respond to deviations from the long-term equilibrium
#' 
#' @param x object of class \code{\link{VECM}} or \code{\link{TVECM}}
#'@param add.legend logical. Whether to add a legend?
#'@param legend.location character. Location of the legend, see \code{\link[graphics]{legend}}
#'@param \ldots arguments passed to the initial \code{plot} call, see \code{\link[graphics]{plot}}
#'@examples
#'data(zeroyld)
#'vec_l1 <- VECM(zeroyldMeta[, c("long.run", "short.run")], lag =1)
#'tvec_l1 <- TVECM(zeroyldMeta[, c("long.run", "short.run")], lag =1, 
#'                 plot  = FALSE, trace = FALSE, th1 = list(exact = -1.263))
#'
#'plot_ECT(vec_l1)
#'plot_ECT(tvec_l1, legend.location = "bottomright")
#'@export
#'@return a plot, and invisibly the underlying data.frame, containing the ECT and the response for each variable

plot_ECT <-  function(x, add.legend = TRUE, legend.location = "topright", ...) {
  
  
  if(!inherits(x, c("VECM", "TVECM"))) stop("Only works for VECM/TVECM  object")
  ## prep data
  ect_df <- get_ect_resp(x)
 
  
## plot now  
  plot(NA, xlim= range(ect_df$ect), 
       ylim=range(c(ect_df[, -1])), xlab = "ECT", ylab = "Response", ...)
  invisible <- lapply(2:ncol(ect_df), function(x) lines(x = ect_df$ect, y = ect_df[,x], col =x -1)) 
  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)
  if(add.legend) legend(legend.location, legend=colnames(ect_df[, -1]), fill=1:(nrow(ect_df)-1))
  invisible(ect_df)
}


#### Utility function: get_ect ####
get_ect <-  function(x) {
  if(inherits(x, "VECM")) {
    res <- as.data.frame(x$model)$ECT[ -seq_len(x$lag+1)]    
  } else if(inherits(x, "TVECM")) {
    res <- x$model.specific$ect[, 1]
  } else {
    stop("Does not work for this input")
  }
  res
}



#### Utility function: get_ect_resp ####
get_ect_resp <-  function(x) {
  
  ## prep data
  ect_df <- data.frame(ect =sort(get_ect(x)))
  co <- coef(x)
  series <- get_series(x)
  nthresh <- get_nthresh.nlVar(x)
  
  
  ## function to find right coef
  if(nthresh==0) {
    co_fun <- function(ect, i) co[i, "ECT"] 
  } else if(nthresh==1) {
    co_fun <- function(ect, i) ifelse(ect <= getTh(x),
                                      co$Bdown[i, "ECT"],
                                      co$Bup[i, "ECT"])
  } else if(nthresh==2) {
    co_fun <- function(ect, i) ifelse(ect <= getTh(x),
                                      co$Bdown[i, "ECT"],
                                      co$Bup[i, "ECT"])
  }
  
  ## apply function
  for(i in seq_len(x$k)) {
    ect_df[, series[i]] <- ect_df$ect * co_fun(ect_df$ect, i)
  }
  ect_df
}

if(FALSE) {
  library(tsDyn)
  vec_l1 <- VECM(zeroyldMeta[, c("long.run", "short.run")], lag =1)
  tvec_l1 <- TVECM(zeroyldMeta[, c("long.run", "short.run")], lag =1, plot  = FALSE)
  head(tsDyn:::get_ect(x=tvec_l1))
  
  head(tsDyn:::get_ect_resp(x=tvec_l1))
  head(tsDyn:::get_ect_resp(x=vec_l1))

  plot_ECT(x=vec_l1)
  plot_ECT(x=tvec_l1, legend.location = "bottomright")
}