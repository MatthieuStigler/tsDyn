
#' Long-term mean of an AR(p) process
#' 
#' Computes the long term mean of an AR process
#' @param object an object of class \code{\link{linear}}, \code{\link{setar}} or \code{\link{lstar}}
#' @param \ldots unused argument
#' @details The function computes the long-term mean of an AR(p) process, or of the corresponding sub-regimes in 
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
  res <- x[1]/(1- sum(x[-1]))
  names(res) <-  "ar_mean" # might be overwritten
  res
}

## check parameters, use fo_const if constant
check_inc <-  function(include, coef, fo_const) {
  if(include == "none") {
    return(c("ar_mean"=0))
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
    names(res) <- paste("ar_mean", names(coefs_by_regime), sep = "_")
    res
  }
  
  check_inc(get_include(object), coef_set, ar_mean_const_byregime)
}

#'@rdname ar_mean
#' @export
ar_mean.lstar <-  function(object, ...) ar_mean.setar(object = object, ...)

#' @export
ar_mean.star <-  function(object, ...) ar_mean.setar(object = object, ...)


if(FALSE) {
  library(tsDyn)
  lin_cst_l1 <-  linear(lh, m = 1, include = "const")
  lin_trd_l1 <-  linear(lh, m = 1, include = "trend")
  lin_bth_l1 <-  linear(lh, m = 1, include = "both")
  lin_non_l1 <-  linear(lh, m = 1, include = "none")
  
  lin_cst_l2 <-  linear(lh, m = 2, include = "const")

  ## collect all
  ar_mean(lin_non_l1)
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


###### Characteristic roots
root_oneReg <- function(coef, regime=c("L", "M", "H", "."), lags, warn_root = TRUE){
  regime <- match.arg(regime)
  
  ## get reelvant coefs
  coefName <- paste("phi", regime, sep = "")
  phi <- coef[grep(coefName, names(coef))]
  vec <- rep(0, max(lags))
  vec[lags] <- phi
  charPol <- c(1, -phi)
  if(any(is.na(charPol))) {
    warning("some NA coefficients?")
    charPol <-  charPol[!is.na(charPol)]
  }
  root <- polyroot(charPol)
  root_mod <- Mod(root)
  
  if(any(root_mod<=1) & warn_root){
    regimeName <- switch(regime, "L"="low", "M"="medium", "H"="high", "."="")
    message <- paste("Possible unit root in the", regimeName, " regime. Roots are:", paste(round(root_mod,4), collapse=" "))
    warning(message, call.=FALSE)
  }  
  return(list(root = root, root_mod = root_mod))
}


#' Characteristic roots of the AR coefficients
#' 
#' Computes the (inverse) characteristic roots of the auto-regressive coefficients. 
#' To be stationary, the values should be outside the unit circle. The function here returns the
#' modulus of the roots. 
#' @param object object of class \code{\link{nlar}}
#' @param \ldots currently unused
#' @details Computes the roots of the polynomial (1, -phi) using function \code{\link{polyroot}}
#' @return a data.frame, with the modulus of the roots. 
#' For models with multiple regimes (setar, lstar, star), one column per regime.
#' @examples 
#' mod.ar <-  linear(lh, m = 5, include = "const")
#' mod.setar <-  setar(lh, m = 5, include = "const")
#' 
#' charac_root(mod.ar)
#' charac_root(mod.setar)
#' @export
charac_root <-  function(object, ...) UseMethod("charac_root")

#'@rdname charac_root
#'@export
charac_root.nlar <- function(object, ...) {

  model <- object$model.specific$nthresh
  if(inherits(object, c("lstar", "star"))) model <- 1
  reg <-  switch(as.character(model), 
                 "0" = ".", 
                 "1" = c("L", "H"), 
                 "2" = c("L", "M", "H"))
  reg_names <- paste("value", ifelse(reg==".", "all", reg), sep = "_")
  
  ## lags
  mod <- object$model.specific
  lags <-  switch(as.character(model), 
                  "0" = list(seq_len(object$str$m)), 
                  "1" = list(mod$ML, mod$MH), 
                  "2" = list(mod$ML, mod$MM, mod$MH))
  if(inherits(object, c("lstar", "star"))) lags <- list(1:object$str$m, 1:object$str$m)
  
  res_li  <- mapply(root_oneReg, regime = as.list(reg), lags = lags, 
                    MoreArgs = list(coef = coef(object), warn_root = FALSE),
                    SIMPLIFY = FALSE)
  res_li2 <- lapply(res_li, function(x) data.frame(root = 1:length(x[[1]]), value = x$root_mod))
  res_df <- Reduce(function(x, y) merge(x, y, by = "root", all = TRUE), res_li2)
  colnames(res_df)[-1] <- reg_names
  res_df
}

#'@export
charac_root.aar <- function(object, ...) stop("Not available")
  
if(FALSE) {
  library(tsDyn)
  lin_cst_l5 <-  linear(lh, m = 5, include = "const")
  set_cst_l5 <-  setar(lh, m = 5, include = "const")
  lst_cst_l5 <-  lstar(lh, m = 5, include = "const", trace = FALSE)
  set2_cst_l5 <-  setar(lh, ML = c(1, 4), mH = 5, mM =1, include = "const", nthresh = 2)
  mod.star <- star(lynx, mTh=c(0,1), control=list(maxit=3000), trace = FALSE)
  mod.aar <- aar(log(lynx), m=3)
  
  charac_root(object = lin_cst_l5)
  charac_root(object = set_cst_l5)
  charac_root(object = set2_cst_l5)
  charac_root(object = lst_cst_l5)
  charac_root(mod.star)
  charac_root(mod.aar)
}