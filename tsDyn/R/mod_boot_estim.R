
mod_boot <- function(x) {
  x_b <- switch(class(x)[1],
                "linear" = linear.boot(x),
                "setar" = setar.boot(x),
                "VAR" = VAR.boot(x),
                "TVAR" = TVAR.boot(x))
  x_b
}

mod_boot_estim <-  function(x){
  x_b <-  mod_boot(x)
  
  switch(class(x)[1],
          "linear" = linear(x_b, m = x$str$m, include = x$include),
          "setar" = setar(x_b, m = x$str$m, include = x$include,
                          nthresh = x$model.specific$nthresh, th = getTh(x)),
          "VAR" = lineVar(x_b, lag = x$lag, include = x$include),
          )
  
}


if(FALSE) {
  library(tsDyn)
  mode_ar <- linear(lh, m = 2)
  mode_setar <- setar(lh, m = 2)
  mode_VAR <- lineVar(barry, lag =1)

  mod_all <- list(mode_ar = mode_ar,
                  mode_setar = mode_setar,
                  mode_VAR = mode_VAR)  
  lapply(mod_all, mod_boot_estim)
  
}