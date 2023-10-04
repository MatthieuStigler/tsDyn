#' @importFrom generics tidy
#' @export
generics::tidy



#' @export
tidy.VAR <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  
  ## extract from summary  
  coef_df <- as.data.frame(summary(x)$coefMat)
  colnames(coef_df) <- c("estimate", "std.error", "statistic", "p.value")
  
  ## add equation and term columns
  coef_df$equation <- rep(eqNames(x), each = x$nparB)
  coef_df$term <- rep(colnames(coef(x)), x$k)
  coef_df <- coef_df[, c("equation", "term", "estimate", "std.error", "statistic", "p.value")]
  
  if(conf.int){
    CI <- confint(x)
    if(!all(rownames(CI) == rownames(coef_df))) warning("Problem in matching names, contact maintainer.")
    coef_df$term
    coef_df$conf.low <- CI[,1]
    coef_df$conf.high <- CI[,2]
  }
  
  rownames(coef_df) <- NULL
  coef_df
}


if(FALSE){
  library(tsDyn)
  # data(barry)
  x <- lineVar(barry, lag=3)
  tidy(x)
  
  environment(tidy.VAR) <- environment(VECM)
  tidy.VAR(x)
}