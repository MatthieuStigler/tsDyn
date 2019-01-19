
### GIRF plot
#'@param x output of girf
#'@param n_simu line plot: which simulation to plot?
#'@param var plot: which variable to plot?
#'@param plot_type density for each \code{n.ahead}, or line, for multipe \code{n_simu}?
#' @rdname GIRF
#' @export
plot.GIRF_df <- function(x, plot_type = c("density", "line"),
                         n.ahead = c(1, 5, 10), 
                         var = unique(x$var)[1], 
                         n_simu = c(1,2), ...) {
  plot_type <-  match.arg(plot_type)
  
  if(plot_type == "density") {
    plot_GIRF_dens(x=x, n.ahead = n.ahead, var = var, ...)
  } else if(plot_type == "line") {
    plot_GIRF_line(x=x, n_simu = n_simu, var = var, ...)     
  }
}



#'@importFrom stats reshape
plot_GIRF_dens <- function(x, n.ahead = c(1, 5, 10), var = unique(x$var)[1],  
                           add_legend = TRUE, ...) {
  df <- x[x$n.ahead %in% n.ahead & x$var %in% var, 
          c("n_simu", "var", "n.ahead", "girf")]
  
  df_w <- reshape(df, idvar = c("n_simu", "var"), timevar = "n.ahead", direction = "wide", v.names = "girf")
  colnames(df_w) <- gsub("girf\\.", "n.ahead: ", colnames(df_w))
  
  dens <- apply(df_w[, -c(1, 2), drop = FALSE], 2, density)
  plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")), ylab = "GIRF")
  invisible <- mapply(lines, dens, col=1:length(dens))
  
  if(add_legend) legend("topright", legend=names(dens), fill=1:length(dens))
  
}


plot_GIRF_line <- function(x, n_simu = 1:5, var = unique(x$var)[1], 
                           add_legend = TRUE, ...) {
  df <- x[x$n_simu %in% n_simu & x$var %in% var, 
          c("n_simu", "var", "n.ahead", "girf")]
  
  df_w <- reshape(df, idvar = c("n.ahead", "var"), timevar = "n_simu", direction = "wide", v.names = "girf")
  colnames(df_w) <- gsub("girf\\.", "n_simu: ", colnames(df_w))
  df_w2 <- df_w[, -c(1,2),drop = FALSE]
  
  
  plot(NA, xlim= range(df$n.ahead), 
       ylim=range(df$girf), ylab = "GIRF")
  invisible <- lapply(1:ncol(df_w2), function(x) lines(x = df_w$n.ahead, y = df_w2[,x], col =x))
  if(add_legend) legend("topright", legend=colnames(df_w2), fill=1:nrow(df_w2))
  
}

## lower level function for unexported irf_1_shock
plot_GIRF_line_low <- function(x, n_simu = 1:2,...) {
  if(!"girf" %in% colnames(x)) x$girf <-  x$sim_1 - x$sim_2
  if(!"n_simu" %in% colnames(x)) {
    x$n_simu <-  1
    n_simu <-  1
  }
  if(!"var" %in% colnames(x)) {
    x$var <-  "x"
  }
  plot_GIRF_line(x, n_simu = n_simu, ...) 
}  


