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

## compare
all.equal(coef(all_lin[["const_l_1"]])[2], ar_1_Mean$ar[,,1], check.attributes = FALSE)
all.equal(coef(all_lin[["none_l_1"]]), ar_1_noMean$ar[,,1], check.attributes = FALSE)

all.equal(coef(all_lin[["const_l_2"]])[-1], ar_2_Mean$ar[,,1], check.attributes = FALSE)
all.equal(coef(all_lin[["none_l_2"]]), ar_2_noMean$ar[,,1], check.attributes = FALSE)

## mean?
