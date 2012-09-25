library(tsDyn)

x <- log10(lynx)

### Estimate models
mod <- list()
mod[["linear"]] <- linear(x, m=2)
mod[["setar"]] <- setar(x, m=2, thDelay=1)
mod[["lstar"]] <- lstar(x, m=2, thDelay=1)
mod[["aar"]] <- aar(x, m=2)


### Extract methods
sapply(mod, AIC)
sapply(mod, BIC)
sapply(mod, mse)
sapply(mod, MAPE)

sapply(mod, coef)
sapply(mod, function(x) head(residuals(x)))


lapply(mod, predict, n.ahead=10)