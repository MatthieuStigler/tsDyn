
### Functions to boot, reestimate, and both
## mod_boot: use x.boot
## mod_refit: recall the function, using all input arguments MINUS the user-specified th, etc
## this function could use the call() approach instead, need to think about it!
## mod_boot_estim: do both

## mod_refit_check: internal check


####### mod_boot #######

mod_boot <- function(x) {
  x_b <- switch(class(x)[1],
                "linear" = linear.boot(x),
                "setar" = setar.boot(x),
                "VAR" = VAR.boot(x),
                "VECM" = VECM.boot(x),
                "TVAR" = TVAR.boot(x),
                "TVECM" = TVECM.boot(x),
                stop("ERROR, not implemented"))
  x_b
}

####### mod_refit #######
mod_refit <-  function(x, data) UseMethod("mod_refit")

#'@export
mod_refit.linear <-  function(x, data) {
  x_mod <-  x$model.specific
  linear(data, 
         m = x$str$m, d = x$str$d, steps = x$str$steps, series = x$str$series,
         include = x$include, type = x$x_mod$type)
}

#'@export
mod_refit.setar <-  function(x, data) {
  x_mod <-  x$model.specific
  setar(data, 
        m = x$str$m, d = x$str$d, steps = x$str$steps, series = x$str$series,
        ML = x_mod$ML, MM = x_mod$MM, MH = x_mod$MH, 
        thDelay = x_mod$thDelay, nthresh = x_mod$nthresh, 
        include = x$include, type = x$model.specific$type,
        trim = x_mod$trim, nested = TRUE, 
        trace = FALSE)
}

#'@export
mod_refit.lstar <-  function(x, data) {
  x_mod <-  x$model.specific
  lstar(data, 
        m = x$str$m, d = x$str$d, steps = x$str$steps, series = x$str$series,
        mL = x_mod$mL, mH = x_mod$mH, 
        mTh = x_mod$mTh,
        include = x_mod$include, 
        trace = FALSE)
}

#'@export
mod_refit.VAR <-  function(x, data) {
  x_mod <-  x$model.specific
  inpt_args <- x$inputArgs
  lineVar(data, 
          lag = x$lag, 
          include = x$include, 
          I = inpt_args$I,
          estim = inpt_args$estim,
          model = class(x)[1],
          LRinclude = inpt_args$LRinclude)
}

#'@export
mod_refit.TVAR <-  function(x, data) {
  x_mod <-  x$model.specific
  inpt_args <- x$inputArgs
  TVAR(data, 
       lag = x$lag, 
       include = x$include, 
       thDelay = x_mod$thDelay,
       model = inpt_args$model,
       commonInter = inpt_args$commonInter,
       nthresh = x_mod$nthresh,
       mTh = inpt_args$mTh,
       dummyToBothRegimes = inpt_args$dummyToBothRegimes,
       trim = x$trim, 
       trace = FALSE)
}

#'@export
mod_refit.TVECM <-  function(x, data) {
  x_mod <-  x$model.specific
  inpt_args <- x$inputArgs
  TVECM(data, 
       lag = x$lag, 
       include = x$include, 
       nthresh = x_mod$nthresh,
       dummyToBothRegimes = inpt_args$dummyToBothRegimes,
       trim = inpt_args$trim, 
       trace = FALSE, plot = FALSE)
}


####### mod_refit_check #######

mod_refit_check <- function(object, check_call = FALSE, keep_optim=FALSE) {
  is_uni <-  inherits(object, "nlar")
  if(is_uni) {
    dat <-  object$str$x
  } else {
    dat <-  object$model[, seq_len(object$k)]
  }
  mod_new <- mod_refit(object, data = dat)
  
  if(!check_call) {
    object$inputArgs$call <- NULL
    mod_new$inputArgs$call <- NULL
  }
  if(!keep_optim){
    if("model.specific" %in% names(object) &&"counts" %in% names(object$model.specific)){
      mod_new$model.specific$counts <- NULL
      object$model.specific$counts <- NULL
    }
    
  }
  
  all.equal(object, mod_new)
}


mod_boot_estim <-  function(x){
  x_b <-  mod_boot(x)
  mod_refit(x, x_b)
}


if(FALSE) {
  library(tsDyn)
  mod_ar <- linear(lh, m = 2)
  mod_setar <- setar(lh, m = 2)
  mod_lstar <- lstar(lh, m = 2, trace = FALSE)
  mod_VAR <- lineVar(barry, lag =1)
  mod_VECM <- VECM(barry, lag =1)
  mod_TVAR <- TVAR(barry, lag =1, trace = FALSE)
  mod_TVECM <- TVECM(barry[, 1:2], lag =1, trace = FALSE, plot = FALSE)
  
  ## test
  mod_refit_check(object = mod_ar)
  mod_refit_check(object = mod_setar)
  mod_refit_check(object = mod_lstar)
  mod_refit_check(object = mod_VAR)
  mod_refit_check(object = mod_VECM)
  mod_refit_check(object = mod_TVAR)
  mod_refit_check(object = mod_TVECM)
  
}