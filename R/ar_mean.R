
#' Long-term mean of an AR(p) process
#' 
#' Computes the long term mean of an AR process
#' @param object an objetc of class \code{\link{linear}}, \code{\link{setar}} or \code{\link{lstar}}
#' @param \ldots unused argment
#' @details The function computes the long-term meanof an AR(p) process, or of the correspondg sub-regimes in 
#' SETAR or LSTAR model. 
#' There are three possible cases:
#' \describe{
#'   \item{No constant nor trend}{The LT mean is 0}
#'   \item{constant}{The LT mean is given by const/(1-sum(AR coefs))}
#'   \item{Trend}{The LT mean is not defined}
#'}
#' @examples 
#' ## estimate a (linear) AR, a SETAR and a LSTAR
#'lin_cst_l1 <-  linear(lh, m = 1, include = "const")
#'set_cst_l1 <-  setar(lh, m = 1, include = "const") 
#'lst_cst_l1 <-  lstar(lh, m = 1, include = "const", trace = FALSE)
#'
#'ar_mean(lin_cst_l1)
#'ar_mean(set_cst_l1)
#'ar_mean(lst_cst_l1)

#' @export
ar_mean <- function(object, ...) UseMethod("ar_mean")


### utiliy function ###

## compute the mean if has constant
ar_mean_const <- function(x) {
  x[1]/(1- sum(x[-1]))
}

## check parameters, use fo_const if constant
check_inc <-  function(include, coef, fo_const) {
  if(include == "none") {
    return(0)
  } else if(include == "const") {
    return(fo_const(coef))
  } else if(include %in% c("trend", "both")) {
    warning("No long term mean for an AR with trend")
    return(NULL)
  }
}

#'@rdname ar_mean
#' @export
ar_mean.linear <-  function(object, ...) {
  check_inc(object$include, coef(object), ar_mean_const)
}

#'@rdname ar_mean
#' @export
ar_mean.setar <-  function(object, ...) {
  coef_set <-  coef(object, hyperCoef = FALSE)
  ar_mean_const_byregime <- function(x) {
    reg_co <- gsub('const\\.|phi|\\.[0-9]+$','', names(x))
    coefs_by_regime <- split(x, reg_co)
    res <- sapply(coefs_by_regime, ar_mean_const)
    names(res) <- names(coefs_by_regime)
    res
  }
  
  check_inc(get_include(object), coef_set, ar_mean_const_byregime)
}

#'@rdname ar_mean
#' @export
ar_mean.lstar <-  function(object, ...) ar_mean.setar(object = object, ...)


if(FALSE) {
  library(tsDyn)
  lin_cst_l1 <-  linear(lh, m = 1, include = "const")
  lin_trd_l1 <-  linear(lh, m = 1, include = "trend")
  lin_bth_l1 <-  linear(lh, m = 1, include = "both")
  lin_non_l1 <-  linear(lh, m = 1, include = "none")
  
  lin_cst_l2 <-  linear(lh, m = 2, include = "const")

  ## collect all
  ar_mean(lin_cst_l1)
  lin_all <-  list(lin_cst_l1 = lin_cst_l1,
                   lin_trd_l1 = lin_trd_l1,
                   lin_bth_l1 = lin_bth_l1,
                   lin_non_l1 = lin_non_l1,
                   lin_cst_l2 = lin_cst_l2)

  suppressWarnings(sapply(lin_all, ar_mean))
  
  ## setar
  set_cst_l1 <-  setar(lh, m = 1, include = "const")
  set_trd_l1 <-  setar(lh, m = 1, include = "trend")
  set_bth_l1 <-  setar(lh, m = 1, include = "both")
  set_non_l1 <-  setar(lh, m = 1, include = "none")

  set_all <-  list(set_cst_l1 = set_cst_l1,
                   set_trd_l1 = set_trd_l1,
                   set_bth_l1 = set_bth_l1,
                   set_non_l1 = set_non_l1)

  suppressWarnings(sapply(set_all, ar_mean))
  
  ## lstar
  lst_cst_l1 <-  lstar(lh, m = 1, include = "const", trace = FALSE)
  lst_trd_l1 <-  lstar(lh, m = 1, include = "trend", trace = FALSE)
  lst_bth_l1 <-  lstar(lh, m = 1, include = "both", trace = FALSE)
  lst_non_l1 <-  lstar(lh, m = 1, include = "none", trace = FALSE)
  
  coef(lst_cst_l1)
  
  lst_all <-  list(lst_cst_l1 = lst_cst_l1,
                   lst_trd_l1 = lst_trd_l1,
                   lst_bth_l1 = lst_bth_l1,
                   lst_non_l1 = lst_non_l1)
  get_include(lst_cst_l1)
  ar_mean(lst_cst_l1)
  sapply(lst_all, ar_mean)
  suppressWarnings(sapply(lst_all, ar_mean))
}