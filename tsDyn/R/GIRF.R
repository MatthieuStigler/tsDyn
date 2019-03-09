
## see also: https://www.econstor.eu/bitstream/10419/44961/1/65618079X.pdf on ideas: https://ideas.repec.org/p/zbw/bubdp1/201103.html
## on tsDyn: https://github.com/MatthieuStigler/tsDyn_GIRF/blob/master/README.md
## https://github.com/angusmoore/tvarGIRF

## terasvirta: https://pure.au.dk/ws/files/54473123/rp13_18.pdf


#' Generalized Impulse response Function (GIRF)
#' 
#' Generates a GIRF for multiple innovations and histories
#' 
#' @param object An object of class \code{\link{linear}}, \code{\link{setar}} or \code{nlVar} (\code{\link{TVAR}}, \code{\link{TVECM}})
#' @param n.ahead The numher of steps ahead to compute
#' @param seed optional, the seed for the random numbers
#' @param \ldots Further arguments passed to specific methods. 
#' @details In a nonlinear model, the Impulse response Function (IRF) is not time-, scale- and sign-invariant as in linear models. 
#' To cope with this, Koop et al (1996) introduced the Generalized Impulse response Function (GIRF):
#' \deqn{IRF(h,\delta,\omega_{t-1})=E[y_{t+h}|\epsilon_{t}=\delta,\epsilon_{t+h}=0,\omega_{t-1}]-E[y_{t+h}|\epsilon_{t}=0,\epsilon_{t+h}=0,\omega_{t-1}]}
#' 
#' It is the difference between two conditional expectations, one containing the shock of interest, the second one averaging it out. The averaging-out is done
#' by comparing against random innovations (unlike the IRF, that compare against innovation 0), 
#' The parameter \code{R} corresponds to the number of times this is done. 
#' 
#' The GIRF as defined here depends on the particular shock, as well as history. 
#' Koop et al (1996) suggest to draw multiple combinations of histories and innovations. This is done with arguments \code{n.hist} and \code{n.shock} 
#' (or, alternatively, provide one of, or both, \code{hist_li} and \code{hist_li} as list of histories and shocks).
#' 
#' The output is a data-frame containing the two average paths for each combinations of shocks and histories. 
#'@return A data-frame, with:
#'\describe{ 
#'\item{n_simu:}{Id for the simulation (total number is n.hist times n.shock)} 
#'\item{hist, shock}{History and shock used in the nth simulation}
#'\item{n.ahead:}{The forecasting horizon. Note the shocks happens at time 0}
#'\item{var:}{The variable (on which the shock happens, corresponds hence to the \code{response} argument in \code{irf})}
#'\item{sim_1, sim_2}{The average (over R times) simulation with the specific shock (sim_1) or with ranodm shocks (sim_2). }
#'\item{girf}{The difference between sim_1 and sim_2}
#'}
#'@author Matthieu Stigler
#'@seealso \code{\link{irf.nlVar}} for the IRF, for linear models, or in case of non-linear models, for each regime. 
#'@examples 
#'
#'## simulate a SETAR for the example. Higher regime more persistent (AR coef 0.5 instead of 0.2)
#'set <- setar.sim(B = c(0.5, 0.2, 0.5, 0.5), lag = 1, Thresh = 0.5, n = 500)
#'set_estim <- setar(set, m = 1)
#'
#'## regime-specific IRF:
#'plot(irf(set_estim, regime = "L", boot = FALSE))
#'plot(irf(set_estim, regime = "H", boot = FALSE))
#'
#'## GIRF
#'girf_out <- GIRF(set_estim)
#'
#'## the GIRF shows a very fast convergence (the shock at n.ahead = 4 is already very close to 0)
#'plot(girf_out, n.ahead = 1:4)
#'## investigate a few specific GIRFS:
#'plot(girf_out, plot_type = "line", n_simu  = 1:5)

#'@export
GIRF <-  function(object, n.ahead, seed = NULL, ...) UseMethod("GIRF")

### This function generates, for a given shock and hist, parallel IRF paths, once  with shock, once without. 
### here, innov is input, is missing, ONE set is drawn
irf_1_shock <-  function(object, shock, hist, n.ahead=10, innov= NULL, shock_both = TRUE,
                         returnStarting = FALSE, 
                         add.regime = FALSE, 
                         seed = NULL) {
  

  ## extract model infos
  lag <- get_lag(object)
  n_start <- ifelse(inherits(object, c("VECM", "TVECM")), lag + 1, lag)
  include <- object$include
  B <- coef(object, hyperCoef = FALSE)
  nthresh <-  get_nthresh(object)
  Thresh <- getTh(object)
  N <- n.ahead + 1
  K <-  get_nVar(object)
  series <-  get_series(object)
  beta <- object$model.specific$coint
  
  ## model
  model <- switch(class(object)[[1]],
                  "ar" = "setar",
                  "linear" = "setar",
                  "setar" = "setar",
                  "VAR" = "TVAR",
                  "TVAR" = "TVAR",
                  "VECM" = "TVECM",
                  "TVECM" = "TVECM",
                  stop("Error model not recognised!"))
  
  if(class(object)[[1]] %in% c("TVAR", "TVECM")) {
    B <-  do.call("cbind", B)
  }
  
  ## sample innov, if not provided
  if(is.null(innov)) {
    res_obj <- as.matrix(residuals(object))
    res_obj <- res_obj[-seq_len(lag),, drop = FALSE]
    if(!is.null(seed)) set.seed(seed)
    index_samp <- sample(seq_len(nrow(res_obj)), N, replace = FALSE)
    innov <-  res_obj <- res_obj[index_samp,, drop = FALSE]
  }
  
  ## hist: as matrix
  hist_M <- as.matrix(hist)
  shock_M <- as.matrix(shock)
  
  ## check args
  if(nrow(hist_M)!= n_start) stop("hist should be of same nrow as lag (+1 if VECM)")
  if(ncol(hist_M)!= K) stop("hist should be of same ncol as number of variables")
  if(nrow(shock_M)!= 1) stop("shock should have only one row")
  if(ncol(shock_M)!= K) stop("shock should be of same ncol as number of variables")
  if(nrow(innov)!=N) stop("innov should be of same length as n.ahead + 1")
  
  ## shocks
  innov_1 <- rbind(shock_M, 
                   innov[-1,, drop = FALSE])
  if(shock_both) innov_1 <- innov + rbind(shock_M, matrix(0, nrow=n.ahead, ncol = K))
  innov_2 <- innov
  
  ## steps 3 and 4 in Koop et al (1996)
  sim_1 <- model.gen(model = model,
                     B = B, lag = lag, include= include,
                     nthresh = nthresh, Thresh = Thresh,
                     beta = beta,
                     starting = hist_M,
                     innov = innov_1, n = N,
                     returnStarting = TRUE, add.regime = add.regime)
  sim_2 <- model.gen(model = model,
                     B = B, lag = lag, include= include,
                     nthresh = nthresh, Thresh = Thresh,
                     beta = beta,
                     starting = hist_M,
                     innov = innov_2, n = N,
                     returnStarting = TRUE, add.regime = add.regime)
  

  # assemble data
  n.aheads_all <- c(-n_start:0, seq_len(n.ahead))
  if(K>1) {
    df <- data.frame(var = rep(series, each = nrow(sim_1)),
                     n.ahead = n.aheads_all,
                     sim_1 = unlist(sim_1),
                     sim_2 = unlist(sim_2))
    
    rownames(df) <-  seq_len(nrow(df))
  } else {
    df <- data.frame(var = series,
                     n.ahead = n.aheads_all,
                     sim_1 = sim_1$res,
                     sim_2 = sim_2$res)  
  }
  
  
  
  if(add.regime) {
    df <- df[, 1:4]  
    colnames(df) <-  c("n.ahead", "sim_1", "regime_1", "sim_2")
  }
  
  if(!returnStarting) df <- subset(df, n.ahead>=0)
  # df$diff <- df$sim_1 - df$sim_2
  df
}

### This function runs irf_1_shock R times, drawing R innov sets. Still given shock and hist 
### Output: average IRF
irf_1_shock_ave <- function(object, shock, hist, R=10, n.ahead=10, innov= NULL, shock_both = TRUE,
                            returnStarting = FALSE, 
                            add.regime = FALSE, 
                            seed = NULL) {
  
  if(!is.null(seed)) set.seed(seed)
  out <- replicate(R, irf_1_shock(object = object, shock = shock, hist = hist, n.ahead = n.ahead, 
                                  innov = innov, shock_both  =shock_both,
                                  returnStarting = returnStarting, add.regime = add.regime), simplify = FALSE)
  out_M <- do.call("rbind", out) 
  # out_M$repli = rep(seq_len(R), each = length(unique(out_M$n.ahead))) 
  
  ## step 5 in Koop et al:
  out_M_means <- aggregate(out_M[, grep("sim|regime", colnames(out_M))],
                           list(n.ahead = out_M$n.ahead, var = out_M$var), mean)
  out_M_means
}



#' @rdname GIRF
#' @param n.hist The number of past histories to consider. Should be high, ideally size of data (minus lags). 
#' @param n.shock The number of actual shocks to consider
#' @param R the number of draws to use for the n.ahead innovations
#' @param hist_li optional, a list of histories (each of same length as lags in the model). 
#' If not provided, \code{n.hist} histories will be randomly drawn for the original series. 
#' @param shock_li optional, a list of innovations. 
#' If not provided, \code{n.shock} shocks will be randomly drawn from the estimated residuals.  
#' @references Koop, G, Pesaran, M. H. & Potter, S. M. (1996) Impulse response analysis in nonlinear multivariate models. Journal of Econometrics, 74, 119-147  
#' @export
## GIRF uses irf_1_shock_ave for a large number of draws of shock and histories
GIRF.setar <-  function(object, n.ahead = 10, seed = NULL, n.hist=20, n.shock=20, R = 10, 
                        hist_li = NULL, shock_li = NULL, ...) {

  ##
  lag <- get_lag(object)
  n_start <- ifelse(inherits(object, c("VECM", "TVECM")), lag + 1, lag) ## VECM requires one more start value as in diff
  x_orig <- get_data_orig(object, as.df = TRUE)
  N <-  nrow(x_orig)
  resids <-  as.data.frame(residuals(object, initVal = FALSE))
  n_used <- nrow(resids)
  nVar <- get_nVar(object)

  
  ## construct hist_li if not provided
  sample_hist <-  function() {
    hist_M <- sample(n_start:N, size = 1, replace = FALSE)
    (hist_M - n_start+ 1) : hist_M
  }
  nrow_length <- function(x) {
    nr <- nrow(x)
    if(is.null(nr)) nr <- length(x)
    nr
  }
    
  if(is.null(hist_li)) {
    if(!is.null(seed)) set.seed(seed)
    # hist_li <- replicate(n.hist, x_orig[sample_hist(),, drop = FALSE], simplify = FALSE)
    samples_hist <- sample(n_start:N, size = n.hist, replace = FALSE)
    hist_li <- lapply(samples_hist, function(i) x_orig[(i - n_start+ 1) : i,, drop = FALSE])
  } else {
    if(!is.list(hist_li)) stop("hist_li should be a list of vectors/matrices")
    if(unique(sapply(hist_li, nrow_length))!=n_start) stop("each element of hist_li should have length lags (+1 if VECM)")
  }
  
  ## construct shock_li if not provided
  if(is.null(shock_li)) {
    if(!is.null(seed)) set.seed(seed)
    samples_shock <- sample(seq_len(n_used), size = n.shock, replace = FALSE)
    shock_li <- lapply(samples_shock, function(i) resids[i,, drop = FALSE])
    # shock_li <- replicate(n.shock, resids[sample(seq_len(n_used), size = 1),, drop = FALSE], simplify = FALSE)
  } else {
    if(!is.list(shock_li)) stop("shock_li should be a list of vectors")
    if(unique(sapply(shock_li, nrow_length))!= 1) stop("each element of shock_li should have length lags")
  }
  
  ## combine each
  M <- expand.grid(hist =hist_li, shock =shock_li)
  n_cases <- nrow(M)
  shock_M <- as.data.frame(do.call("rbind", M$shock))
  colnames(shock_M) <- paste("shock_var", seq_len(nVar), sep = "")
  hist_M <- as.data.frame(do.call("rbind", lapply(M$hist, function(x) c(as.matrix(x)))))
  colnames(hist_M) <- paste("hist_x", rep(seq_len(nVar), each = n_start),
                            "_l", rep(seq_len(n_start), times = nVar), sep = "")
  rep_info <- cbind(n_simu = seq_len(nrow(M)),
                    hist_M, shock_M)
    
  ## run irf_1_shock_ave for each combo
  sims <- lapply(1:nrow(M), function(i) irf_1_shock_ave(object = object, 
                                                shock = M$shock[[i]], hist = M$hist[[i]],
                                                n.ahead = n.ahead, 
                                                R = R, seed = seed, ...))
  
  ## extend the data frame
  sims_df <- do.call("rbind", sims)
  n.ahead_here <- length(unique(head(sims_df$n.ahead, 2*(n.ahead*lag)))) # just in case was called with returnStarting
  sims_df$n_simu <- rep(seq_len(n_cases), each = n.ahead_here * nVar)
  sims_df$girf <- with(sims_df, sim_1 - sim_2)
  sims_df2 <-  merge(rep_info, sims_df, by = "n_simu", sort = FALSE)
  
  ##
  # if("regime_1" %in% colnames(sims_df)) {
  #   cols <- c("n_rep", "last_hist", "shock", "n.ahead", "sim_1", "sim_2", "regime_1", "girf")
  # } else {
  #   cols <- c("n_rep", "last_hist", "shock", "n.ahead", "sim_1", "sim_2", "girf")
  # }
  # sims_df[, cols]
  class(sims_df2) <- c("GIRF_df", "data.frame")
  sims_df2
}

#' @rdname GIRF
#' @export
GIRF.linear <- function(object, n.ahead = 10, seed = NULL, n.hist=20, n.shock=20, R = 10, 
                        hist_li = NULL, shock_li = NULL, ...) {
  
  GIRF.setar(object, n.ahead = n.ahead, seed = seed, n.hist=n.hist, n.shock=n.shock, R = R, 
                          hist_li = hist_li, shock_li = shock_li, ...)
}

#' @rdname GIRF
#' @export
GIRF.nlVar <- function(object, n.ahead = 10, seed = NULL, n.hist=20, n.shock=20, R = 10, 
                        hist_li = NULL, shock_li = NULL, ...) {
  
  GIRF.setar(object, n.ahead = n.ahead, seed = seed, n.hist=n.hist, n.shock=n.shock, R = R, 
             hist_li = hist_li, shock_li = shock_li, ...)
}


if(FALSE) {
  library(tsDyn)
  library(tidyverse)
  
  library(devtools)
  load_all()
  
  set_l2_const <- setar(lh, m = 2, include = "const")
  linear_l2_none <- linear(lh, m = 2, include = "none")
  linear_l2_const <- linear(lh, m = 2, include = "const")
  
  environment(irf_1_shock) <-  environment(setar)
  environment(GIRF.setar) <-  environment(setar)
  irf_1_shock <-  tsDyn:::irf_1_shock
  irf_1_shock_ave <- tsDyn:::irf_1_shock_ave
  
  ##  irf_1_shock
  innov_obs <- residuals(linear_l2_none)
  irf_1_shock(object = linear_l2_none, 
              shock = innov_obs[5],
              hist = linear_l2_none$str$x[c(8, 9)],
              innov = sample(innov_obs[-c(1, 2)], size =11))
  
  ##  irf_1_shock, returnStarting
  irf_1_shock(object = linear_l2_none, 
              shock = innov_obs[5],
              hist = linear_l2_none$str$x[c(8, 9)],
              innov = sample(innov_obs[-c(1, 2)], size =11),
              returnStarting = TRUE, add.regime = TRUE)
  
  ## irf_1_shock_ave
  irf_1_shock_ave(object = linear_l2_none, 
                  shock = innov_obs[5],
                  hist = linear_l2_none$str$x[c(8, 9)]) %>% 
    mutate(diff = sim_1 - sim_2)

  ## irf_1_shock_ave
  irf_1_shock_ave(object = linear_l2_none, 
                  shock = innov_obs[5],
                  hist = linear_l2_none$str$x[c(8, 9)],
                  returnStarting = TRUE, add.regime = TRUE) %>% 
    mutate(diff = sim_1 - sim_2)
  
  ##
  girf_lin_l2 <- GIRF(linear_l2_none, add.regime  =TRUE, returnStarting = FALSE, R = 5) %>%  as_tibble
  filter(girf_lin_l2, n_rep == 1)
  
  girf_lin_l2_df <- as.data.frame(girf_lin_l2)
  
  ## simple irf
  irf_linear_l2 <- irf(linear_l2_none)
  plot(irf_linear_l2)
  plot(irf(linear_l2_none, n.ahead = 50))
  
  ## one specific profile: 8 shocks
  girf_lin_l2 %>% 
    filter(n_rep %in% 1:8) %>% 
    qplot(x = n.ahead, y = girf, colour = factor(n_rep), geom = "line", data = .) +
    geom_hline(yintercept = 0)
  
  ## compare with irf
  girf_lin_l2 %>% 
    filter(n_rep ==1) %>% 
    mutate(girf2 = girf/girf[1]) %>% 
    mutate(irf  = tsDyn:::irf_1.linear(linear_l2_none, n.ahead =11),
           diff_irf_girf = irf - girf2)
  
  ## densities
  girf_lin_l2 %>% 
    filter(n.ahead %in% c(0, 1, 10)) %>% 
    qplot(x = girf, colour = factor(n.ahead), geom = "density", data = .) +
    facet_grid(n.ahead~.)
  
  ## why all same??
  plot(density(linear_l2_none$str$x))
  plot(density(residuals(linear_l2_none, initVal = FALSE)))
  plot(density(residuals(linear_l2_const, initVal = FALSE)))
  girf_lin_l2 %>% 
    arrange(n.ahead) %>% 
    group_by(n.ahead) %>% 
    summarise(mean = mean(girf))
  
  
  ### SETAR ###
  
  ## simple irf
  set_l2_c_irf_1_L <- tsDyn:::irf_1.setar(set_l2_const, regime = "L")
  set_l2_c_irf_1_H <- tsDyn:::irf_1.setar(set_l2_const, regime = "H")
  
  set_l2_c_irf_1_df <- tibble(n.ahead = set_l2_c_irf_1_L$n.ahead,
                              reg_L = set_l2_c_irf_1_L$x, 
                              reg_H = set_l2_c_irf_1_H$x)
  
  set_l2_c_irf_1_df_l <-  set_l2_c_irf_1_df %>% 
    gather(regime, irf, starts_with("reg"))
  qplot(x = n.ahead, y= irf, geom="line", colour = regime, data = set_l2_c_irf_1_df_l)
  irf_1_shock(set_l2_const, shock = 1, hist = c(0, 0))
  
  ## girf
  set_girf <- GIRF(object=set_l2_const,
                   hist_li = list(c(1.6, 1.6), c(2.05, 2.05), c(3.2, 3.3)),
                   shock_li = list(0.01), R = 100) %>% as_tibble
  
  set_girf %>% 
    filter(n_rep %in% 1:8) %>% 
    qplot(x = n.ahead, y = girf, colour = factor(last_hist), geom = "line", data = .) +
    geom_hline(yintercept = 0) +
    ggtitle("last_hist: 25 (upper regime)") +
    scale_x_continuous(breaks = seq(0, 10, by = 2))
  
  ### SETAR
  set_girf_full <- GIRF(object=set_l2_const, R = 5, add.regime = FALSE,
                        returnStarting = TRUE) %>%  as_tibble
  set_girf_1series <- filter(set_girf_full, n_rep ==1) %>% 
    mutate(reg_recalc = regime(set_l2_const, serie = sim_1))
  set_girf_1series
  
  set_girf_full %>% 
    count(regime_1)
  
  
  ## re classify regime
  set_girf_full
  
  
  set_girf_full %>% 
    filter(n.ahead %in% 1:2) %>% 
    arrange(last_hist)

  
  set_girf_full %>% 
    filter(n.ahead ==0) %>% 
    ggplot(aes(x = girf, colour = factor(n.ahead)))+
    geom_density()
  
  set_girf_full %>% 
    filter(n.ahead %in% c(0, 3, 4, 6, 8, 10)) %>% 
    ggplot(aes(x = girf, colour = factor(n.ahead)))+
    geom_density()
}
