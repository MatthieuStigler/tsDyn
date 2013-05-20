\name{accuracy_stat}
\alias{accuracy_stat}
\alias{accuracy_stat.default}
\alias{accuracy_stat.pred_roll}
\title{Forecasting accuracy measures.}
\usage{
  accuracy_stat(object, ...)

  \method{accuracy_stat}{default} (object, true, ...)

  \method{accuracy_stat}{pred_roll} (object, w, ...)
}
\arguments{
  \item{object}{A data-frame, matrix, or object of class
  \code{pred_roll}}

  \item{true}{If \code{object} is just a matrix or
  data-frame, true values to be compared to should be
  supplied}

  \item{w}{Optional. For objects of class \code{pred_roll}
  containing multiple variables, user can specify the way
  to aggregate the specific x-step-ahead into the
  \sQuote{all} category}

  \item{\dots}{Not used currently.}
}
\value{
  A data-frame containing the forecasting accuracy
  measures.
}
\description{
  Compute forecasting accuracies. This is very similar ot
  the \code{\link[forecast]{accuracy}} method form
  \pkg{forecast}.
}
\details{
  The function works either for a simple data.frame or for
  objects \code{pred_roll}. For simple data.frames, the
  argument \code{true}, i.e. a data frame containing the
  true values, has to be provided. For \code{pred_roll}
  objects, the true values are contained in the object, so
  no need (nor possibility) to provide the true values.
}
\examples{
## univariate:
mod_ar <- linear(lynx[1:100], m=1)
mod_ar_pred <- predict_rolling(mod_ar, newdata=lynx[101:114])
accuracy_stat(object=mod_ar_pred$pred, true=mod_ar_pred$true)

## multivariate
data(barry)
mod_var <- lineVar(barry, lag=1)

mod_var_pred <-predict_rolling(object=mod_var, nroll=10, n.ahead=1:3)
accuracy_stat(object=mod_var_pred)
accuracy_stat(object=mod_var_pred, w=c(0.7, 0.2, 0.1))
}
\author{
  Matthieu Stigler
}
\keyword{ts}

