
## see also: https://www.econstor.eu/bitstream/10419/44961/1/65618079X.pdf on ideas: https://ideas.repec.org/p/zbw/bubdp1/201103.html
## on tsDyn: https://github.com/MatthieuStigler/tsDyn_GIRF/blob/master/README.md
## https://github.com/angusmoore/tvarGIRF

## terasvirta: https://pure.au.dk/ws/files/54473123/rp13_18.pdf


#' Generalized Impulse response Function (GIRF)
#' 
#' Generates a GIRF for multiple innovations and histories
#' 
#' @param object An object of class \code{\link{linear}} or \code{\link{setar}}
#' @param n.ahead The numher of steps ahead to compute
#' @param seed optional, the seed for the random numbers
#' @param \ldots Further argumetns passed to specific methods. 
#' @export
GIRF <-  function(object, n.ahead, seed = NULL, ...) UseMethod("GIRF")

### This function generates, for a given shock and hist, parallel IRF paths, once  with shock, once without. 
### here, innov is input, is missing, ONE set is drawn
irf_1_shock <-  function(object, shock, hist, n.ahead=10, innov= NULL, shock_both = TRUE,
                         returnStarting = FALSE, 
                         add.regime = FALSE) {
  

  ## extract model infos
  lag <- object$str$m
  include <- object$include
  B <- coef(object, hyperCoef = FALSE)
  nthresh <-  object$model.specific$nthresh
  if(nthresh >0)  Thresh <- getTh(object)
  N <- n.ahead + 1
  
  ## 
  if(is.null(innov)) {
    res_obj <- residuals(object)[-seq_len(lag)]
    innov <-  sample(res_obj, N, replace = FALSE)
  }
  
  ## check args
  if(length(hist)!= lag) stop("hist should be of same length as lag")
  if(length(innov)!=N) stop("innov should be of same length as n.ahead + 1")
  
  ## shocks
  innov_1 <- c(shock, innov[-1])
  if(shock_both) innov_1 <- innov + c(shock, rep(0, n.ahead))
  innov_2 <- innov
  
  ## steps 3 and 4 in Koop et al (1996)
  sim_1 <- setar.gen(B = B, lag = lag, include= include,
                     nthresh = nthresh, Thresh = Thresh,
                     starting = hist,
                     innov = innov_1, n = N,
                     returnStarting = TRUE, add.regime = add.regime)
  sim_2 <- setar.gen(B = B, lag = lag, include= include,
                     nthresh = nthresh, Thresh = Thresh,
                     starting = hist,
                     innov = innov_2, n = N,
                     returnStarting = TRUE, add.regime = add.regime)
  

  # assemble data
  df <- data.frame(n.ahead = c(-lag:0, seq_len(n.ahead)),
                   sim_1 = sim_1,
                   sim_2 = sim_2)
  if(add.regime) {
    df <- df[, 1:4]  
    colnames(df) <-  c("n.ahead", "sim_1", "regime_1", "sim_2")
  }
  if(!returnStarting) df <- df[-seq_len(lag),]
  # df$diff <- df$sim_1 - df$sim_2
  df
}

### This function runs irf_1_shock R times, drawing R innov sets. Still given shock and hist 
### Output: average IRF
irf_1_shock_ave <- function(object, shock, hist, R=10, n.ahead=10, innov= NULL, shock_both = TRUE,
                            returnStarting = FALSE, 
                            add.regime = FALSE) {
  out <- replicate(R, irf_1_shock(object = object, shock = shock, hist = hist, n.ahead = n.ahead, 
                                  innov = innov, shock_both  =shock_both,
                                  returnStarting = returnStarting, add.regime = add.regime), simplify = FALSE)
  out_M <- do.call("rbind", out) 
  out_M$repli = rep(seq_len(R), each = length(unique(out_M$n.ahead))) 
  
  ## step 5 in Koop et al:
  out_M_means <- aggregate(out_M[, grep("sim|regime", colnames(out_M))],
                           list(n.ahead=  out_M$n.ahead), mean)
  out_M_means
}



#' @rdname GIRF
#' @param n.hist The number of past histories ot consider
#' @param n.shock The number of actual shocks ot consider
#' @param R the number of draws to use for the n.ahead innovations
#' @param hist_li optional, a list of histories (each of same length as lags in the model)
#' @param shock_li optional, a list of innovations
#' @export
## GIRF uses irf_1_shock_ave for a large number of draws of shock and histories
GIRF.setar <-  function(object, n.ahead = 10, seed = NULL, n.hist=20, n.shock=20, R = 10, 
                        hist_li = NULL, shock_li = NULL, ...) {

  ##
  lag <- object$str$m
  x_orig <- object$str$x
  N <-  length(x_orig)
  resids <-  residuals(object, initVal = FALSE)
  
  ## construct list of hist and shocks
  
  ## construct hist_li if not provided
  sample_hist <-  function() {
    hist_M <- sample(lag:N, size = 1, replace = FALSE)
    (hist_M - lag+ 1) : hist_M
  }
  if(is.null(hist_li)) {
    hist_li <- replicate(n.hist, sample_hist(), simplify = FALSE)  
  } else {
    if(!is.list(hist_li)) stop("hist_li should be a list of vectors")
    if(unique(sapply(hist_li, length))!=lag) stop("each element of hist_li should have length lags")
  }
  
  ## construct shock_li if not provided
  if(is.null(shock_li)) {
    shock_li <- replicate(n.shock, sample(resids, size = 1), simplify = FALSE)
  } else {
    if(!is.list(shock_li)) stop("shock_li should be a list of vectors")
    if(unique(sapply(shock_li, length))!= 1) stop("each element of shock_li should have length lags")
  }
  
  ## combine each
  M <- expand.grid(hist =hist_li, shock =shock_li)
  
  ## run irf_1_shock_ave for each combo
  sims <- lapply(1:nrow(M), function(x) irf_1_shock_ave(object = object, 
                                                shock = M$shock[[x]], hist = M$hist[[x]],
                                                n.ahead = n.ahead, 
                                                R = R, ...))
  
  ## extend the data frame
  sims_df <- do.call("rbind", sims)
  n.ahead_here <- length(unique(head(sims_df$n.ahead, 2*(n.ahead*lag)))) # just in case was called with returnStarting
  sims_df$n_rep <- rep(1:nrow(M), each = n.ahead_here)
  sims_df$shock <-  rep(unlist(M$shock), each = n.ahead_here)
  sims_df$last_hist <-  rep(sapply(M$hist, tail, 1), each = n.ahead_here)
  sims_df$girf <- with(sims_df, sim_1 - sim_2)
  
  ##
  if("regime_1" %in% colnames(sims_df)) {
    cols <- c("n_rep", "last_hist", "shock", "n.ahead", "sim_1", "sim_2", "regime_1", "girf")
  } else {
    cols <- c("n_rep", "last_hist", "shock", "n.ahead", "sim_1", "sim_2", "girf")
  }
  sims_df[, cols]
}

#' @rdname GIRF
#' @export
GIRF.linear <- function(object, n.ahead = 10, seed = NULL, n.hist=20, n.shock=20, R = 10, 
                        hist_li = NULL, shock_li = NULL, ...) {
  
  GIRF.setar(object, n.ahead = n.ahead, seed = seed, n.hist=n.hist, n.shock=n.shock, R = R, 
                          hist_li = hist_li, shock_li = shock_li, ...)
}

if(FALSE) {
  library(tsDyn)
  library(tidyverse)
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
  
  set_l2_c_irf_1_df <- tibble(n.ahead = 1:10,
                                  reg_L = set_l2_c_irf_1_L, 
                                  reg_H = set_l2_c_irf_1_H)
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
  set_girf_full <- GIRF(object=set_l2_const, R = 5, add.regime = TRUE,
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
