library(tsDyn)

### linear
grid <-  expand.grid(include = c( "const", "trend","none", "both"),
                     lag = 1:3)
all_lin <- mapply(linear, include = as.character(grid$include),  m = grid$lag, MoreArgs = list(x = lh),
                  SIMPLIFY = FALSE)
names(all_lin) <-  paste(grid$include, "l", grid$lag, sep="_")
lapply(all_lin, coef)

## compare with ar?
ar_1_noMean <- ar(lh, order.max =1, demean = FALSE, method = "ols")
ar_1_Mean <- ar(lh, order.max =1, demean = TRUE, method = "ols")

ar_2_noMean <- ar(lh, aic = FALSE, order.max =2, demean = FALSE, method = "ols")
ar_2_Mean <- ar(lh, aic = FALSE, order.max =2, demean = TRUE, method = "ols")

## compare coefs
all.equal(coef(all_lin[["const_l_1"]])[2], ar_1_Mean$ar[,,1], check.attributes = FALSE)
all.equal(coef(all_lin[["none_l_1"]]), ar_1_noMean$ar[,,1], check.attributes = FALSE)

all.equal(coef(all_lin[["const_l_2"]])[-1], ar_2_Mean$ar[,,1], check.attributes = FALSE)
all.equal(coef(all_lin[["none_l_2"]]), ar_2_noMean$ar[,,1], check.attributes = FALSE)

## compare means: intercept in ar is (1-phi)(diff means)
## diff means is difference between estimated mean, and c/(1- phi)
linear_1_Mean <-  all_lin[["const_l_1"]] 
comp_means <- function(ar, linear) {
  ar_empimean <- ar$x.mean  
  ar_int <- ar$x.intercept  
  lin_coefs <- coef(linear)
  lin_LR_mean <-  ar_mean(linear)
  diff_means <- lin_LR_mean - ar_empimean
  res <- c(ar_int = ar_int, diff = (1- sum(lin_coefs[-1])) * diff_means)
  list(res, check = all.equal(res[1], res[2], check.attributes = FALSE))
}
comp_means(ar = ar_1_Mean, linear = linear_1_Mean)$check
comp_means(ar = ar_2_Mean, linear = all_lin[["const_l_2"]] )$check


## predict?
all_lin_const <- all_lin[grep("const", names(all_lin))]
all_lin_pred <- lapply(all_lin_const, predict, n.ahead = 5)

ar_1_Mean_pred <- predict(ar_1_Mean, n.ahead = 5)$pred

all.equal(ar_1_Mean_pred, all_lin_pred[["const_l_1"]], check.attributes = TRUE)
