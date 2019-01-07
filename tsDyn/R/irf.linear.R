irf_1 <- function(x, n.ahead = 10, cumulative = FALSE, ...) UseMethod("irf_1")

irf_1.ar <-  function(x, n.ahead=10, cumulative=FALSE, ...) {
  coefs <- as.numeric(x$ar)
  empty_series <- c(1, rep(0, n.ahead-1))
  res <- as.numeric(stats::filter(empty_series, coefs, method = "recursive"))
  if(cumulative) res <-  cumsum(res)
  res
}

irf_1.linear <-  function(x, n.ahead=10, cumulative=FALSE, ...) {
  coefs <- coef(x)
  if(any(grepl("const|trend", names(coefs)))) coefs <-  coefs[-grep("const|trend", names(coefs))]
  
  empty_series <- c(1, rep(0, n.ahead-1))
  
  res <- as.numeric(stats::filter(empty_series, coefs, method = "recursive"))
  if(cumulative) res <-  cumsum(res)
  res
}

irf_1.setar <-  function(x, n.ahead=10, cumulative=FALSE, regime = c("L", "M", "H"), ...) {
  regime <-  match.arg(regime)
  coefs <- coef(x, regime = regime, hyperCoef = FALSE)
  if(any(grepl("const|trend", names(coefs)))) coefs <-  coefs[-grep("const|trend", names(coefs))]
  
  empty_series <- c(1, rep(0, n.ahead-1))
  
  res <- as.numeric(stats::filter(empty_series, coefs, method = "recursive"))
  if(cumulative) res <-  cumsum(res)
  res
}


## 
irf_any <-  function(x, n.ahead = 10, cumulative = FALSE, 
                     boot = TRUE, ci = 0.95, runs = 100,
                     regime = c("L", "M", "H"),
                     seed = NULL, 
                     ...) {
  regime <-  match.arg(regime)
  irf_orig <- irf_1(x=x , n.ahead = n.ahead, cumulative = cumulative, regime = regime)
  if(runs ==0) boot <-  FALSE
  
  ## 
  if(boot) {
    if(!is.null(seed)) set.seed(seed) 
    boot_estim <- function(x) {
      x_b <- setar.boot(x)
      if(inherits(x, "linear")) {
        mod_b <- linear(x_b, m = x$str$m, include = x$include)
      } else if(inherits(x, "setar")) {
        mod_b <- setar(x_b, m = x$str$m, include = x$include,
                       nthresh = x$model.specific$nthresh, th = getTh(x))
      }
      
      irf_b <- irf_1(x=mod_b , n.ahead = n.ahead, cumulative = cumulative, regime = regime)
      irf_b
    }
    # boot_estim(x)
    IRF_b <- t(replicate(runs, boot_estim(x)))
    quants <- t(apply(IRF_b, 2, quantile, probs = c(1-ci, ci)))
  }
  
  ## results
  res <- list()
  res$irf$x <- as.matrix(irf_orig)
  colnames(res$irf$x) <-  "x"
  res$model <-  "varest"
  res$runs <- runs
  res$ci <- ci
  if(boot) {
    res$Lower$x <-  matrix(quants[, 1], dimnames = list(NULL, "x"))
    res$Upper$x <-  matrix(quants[, 2], dimnames = list(NULL, "x"))
  }
  res$response <- "x"
  res$ortho <-  FALSE
  res$impulse <- "x"
  res$boot <- boot
  res$cumulative <- cumulative
  
  class(res) <- "varirf"
  res
}

#' @rdname irf.nlVar
#' @export
irf.linear <-  function(x, impulse=NULL, response=NULL, n.ahead=10, ortho=TRUE, cumulative=FALSE, 
                        boot=TRUE, ci=0.95, runs=100, seed=NULL, ...) {
  if(!is.null(impulse) | !is.null(response) | !ortho) stop("Arguments used only for multivariate models")
  irf_any(x=x, n.ahead = n.ahead, cumulative = cumulative, 
           boot = boot, ci = ci, runs = runs, seed = seed, ...)
}

#' @rdname irf.nlVar
#' @param regime For a setar model, which regime (L, M or H) to produce IRF for?
#' @export
irf.setar <-  function(x, impulse=NULL, response=NULL, n.ahead=10, ortho=TRUE, cumulative=FALSE, 
                        boot=TRUE, ci=0.95, runs=100, seed=NULL, regime = c("L", "M", "H"),...) {
  if(!is.null(impulse) | !is.null(response) | !ortho) stop("Arguments used only for multivariate models")
  irf_any(x=x, n.ahead = n.ahead, cumulative = cumulative, 
          boot = boot, ci = ci, runs = runs, regime = regime, seed = seed, ...)
}

#' @rdname irf.nlVar
#' @export
irf.ar <-  function(x, impulse=NULL, response=NULL, n.ahead=10, ortho=TRUE, cumulative=FALSE, 
                    boot=TRUE, ci=0.95, runs=100, seed=NULL, ...) {
  if(!is.null(impulse) | !is.null(response) | !ortho) stop("Arguments used only for multivariate models")
  irf_any(x=x, n.ahead = n.ahead, cumulative = cumulative, 
          boot = boot, ci = ci, runs = runs, seed = seed,  ...)
}

## 
irf_1_sim <-  function(x, n.ahead = 10, innov, ...) {
  
  lag <- x$str$m
  include <- x$include
  B <- coef(x)
  if(any(grepl("const|trend", names(B)))) B[grep("const|trend", names(B))] <-  0
  
  
  start_vals <-  rep(0, lag)
  # shock <- switch(include, "const" = 1 - coef(x)[1],
  #                 "none" = 1)
  if(missing(innov)) innov <-  c(1, rep(0, n.ahead+lag))
  res <- setar.gen(B=B, n=n.ahead, lag=lag, include = include,  
                   nthresh = 0,
                   returnStarting = FALSE,
                   starting=start_vals, innov=innov, ...)
  
  res
}





## example from Sims: https://www3.nd.edu/~esims1/arp_companion.pdf
if(FALSE) {
  coefs <- c(0.8, 0.6, -.5)
  out <- filter(c(1, 0, 0, 0, 0), coefs, method = "recursive")
  out
  sum(out [1:3]* rev(coefs ))
  
}

## TEST
if(FALSE) {
  library(tsDyn)
  ar_2_noMean <- ar(lh, aic = FALSE, order.max = 2, method = "ols", demean = FALSE)  
  ar_2_Mean <- ar(lh, aic = FALSE, order.max =2, demean = TRUE, method = "ols")
  
  linear_l2_none <- linear(lh, m = 2, include = "none")
  linear_l2_const <- linear(lh, m = 2, include = "const")
  setar_l2_const <- setar(lh, m = 2, include = "const")
  all.equal(coef(linear_l2_none), ar_2_noMean$ar[,,1], check.attributes = FALSE)
  all.equal(coef(linear_l2_const)[-1], ar_2_Mean$ar[,,1], check.attributes = FALSE)
  
  ## compar emeans?
  coef(linear_l2_const)
  ar_2_Mean$ar[,,1]
  ar_2_Mean$x.mean
  ar_2_Mean$x.intercept
  
  coefs_ar <- coef(linear_l2_none)
  st_0 <-  1
  st_1 <- coefs_ar[1] * st_0
  st_2 <- coefs_ar[1] * st_1 +  coefs_ar[2] * st_0
  st_3 <- coefs_ar[1] * st_2 +  coefs_ar[2] * st_1
  c(1, st_1, st_2, st_3)
  
  irf_1(ar_2_noMean)
  irf_1(linear_l2_none)
  all.equal(irf_1(x=ar_2_noMean), irf_1(linear_l2_none))
  all.equal(irf_1(x=linear_l2_none),   irf_1_sim(linear_l2_none))
  
  
  ## irf constant
  co_const <- coef(linear_l2_const)

  st_cst_0 <- 1
  st_cst_1 <- co_const[2] * st_cst_0
  st_cst_2 <- co_const[2] * st_cst_1 +  co_const[3] * st_cst_0
  st_cst_3 <- co_const[2] * st_cst_2 +  co_const[3] * st_cst_1
  res_manual <- c(st_cst_0, st_cst_1, st_cst_2, st_cst_3)

  all.equal(res_manual, irf_1_sim(x=linear_l2_const, n.ahead = 4), check.attributes = FALSE)
    
  
  ## irf shock
  # library(tidyverse)
  linear_l2_const_sh <- irf_1_shock(x= linear_l2_const) %>% 
    mutate(n_row = 1:n())
  plot(1:48, linear_l2_const_sh$sim_1, type = "l")
  lines(1:48, linear_l2_const_sh$sim_2, lty = 2, col = 2)
  
  linear_l2_const_sh %>% 
    filter(between(n_row, 11, 20)) %>% 
    qplot(x = n_row, y = diff, data = ., geom = "line")
  
  irf_dat <- data_frame(n_row =11 :(11+10-1), irf = irf_1_sim(x=linear_l2_const))
  qplot(x = n_row, y = irf, data = irf_dat, geom = "line")
  
  ## irf filter
  irf_filt <- irf.linear(x=linear_l2_const)
  all.equal(irf_filt$irf$x, linear_l2_const_sh$diff[11:20])
  
  data_frame(irf_filter = irf.linear(x=linear_l2_const),
             n_row = 1:10) %>% 
    qplot(x = n_row, y = irf_filter, data = ., geom = "line")
  
  
  ### IRF full
  irf_full <-  irf(linear_l2_const)
  
  irf_full$boot <-  FALSE
  irf_full$model <- "varest"
  # fixInNamespace(plot.varirf,  "vars")
  plot(irf_full)
  a(irf_full)
  
  ## setar
  irf_1_L <- irf_1(setar_l2_const, regime = "L")
  irf_1_H <- irf_1(x=setar_l2_const, regime = "H")
  plot(1:10, irf_1_L, type = "l")
  lines(1:10, irf_1_H, lty = 2)
  irf_set_regL <- irf(setar_l2_const, regime = "L", ci = 0.8)
  irf_set_regH <- irf(setar_l2_const, regime = "H", ci = 0.8)
  plot(irf_set_regL, ylim = c(0, 1.2))
  plot(irf_set_regH)
}


