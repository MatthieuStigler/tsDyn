#' @importFrom tseriesChaos embedd mutual
#' @importFrom tseries na.remove
#' @importFrom foreach foreach "%dopar%" "%:%"
#' @importFrom mnormt rmnorm
#' @importFrom graphics abline axis contour curve image layout legend lines lines.default par
#' @importFrom graphics persp plot points segments title
#' @importFrom methods new hasArg
#' @importFrom stats AIC acf aggregate anova ar arima.sim as.formula as.ts coef coefficients
#' @importFrom stats confint confint.default dchisq deltat density deviance embed frequency is.ts lm lm.fit
#' @importFrom stats logLik na.omit optim optimHess pacf pf pgamma plogis plot.ts pnorm predict
#' @importFrom stats princomp printCoefmat pt qgamma quantile residuals rnorm sd symnum time ts
#' @importFrom stats ts.plot tsp "tsp<-" update var vcov weighted.mean window
#' @importFrom utils head tail tail.matrix toLatex

#' @useDynLib tsDyn, .registration = TRUE
NULL