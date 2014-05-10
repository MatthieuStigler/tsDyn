\name{nlar}
\alias{nlar}
\title{Non-linear time series model, base class definition}
\usage{
nlar(str, coefficients, fitted.values, residuals, k, model,
  model.specific = NULL, ...)
}
\arguments{
  \item{str}{a \code{nlar.struct} object, i.e. the result
  of a call to \code{\link{nlar.struct}}}

  \item{coefficients,fitted.values,residuals,k,model,model.specific}{internal
  structure}

  \item{\dots}{further model specific fields}
}
\value{
An object of class \code{nlar}. \link{nlar-methods} for a
list of available methods.
}
\description{
Generic non-linear autogregressive model class constructor.
}
\details{
Constructor for the generic \code{nlar} model class. On a
fitted object you can call some generic methods. For a list
of them, see \code{\link{nlar-methods}}.

An object of the \code{nlar} class is a list of (at least)
components: \describe{ \item{str}{
\code{\link{nlar.struct}} object, encapsulating general
infos such as time series length, embedding parameters,
forecasting steps, model design matrix }
\item{coefficients}{ a named vector of model
estimated/fixed coefficients } \item{k}{ total number of
estimated coefficients } \item{fitted.values}{ model fitted
values } \item{residuals}{ model residuals } \item{model}{
data frame containing the variables used }
\item{model.specific}{ (optional) model specific additional
infos} }

A \code{\link{nlar}} object normally should also have a
model-specific subclass (i.e., \code{nlar} is a virtual
class).

Each subclass should define at least a \code{print} and,
hopefully, a \code{oneStep} method, which is used by
\code{\link{predict.nlar}} to iteratively extend ahead the
time series.
}
\author{
Antonio, Fabio Di Narzo
}
\references{
Non-linear time series models in empirical finance, Philip
Hans Franses and Dick van Dijk, Cambridge: Cambridge
University Press (2000).

Non-Linear Time Series: A Dynamical Systems Approach, Tong,
H., Oxford: Oxford University Press (1990).
}
\seealso{
\code{\link{availableModels}} for currently available
built-in models.  \link{nlar-methods} for available
\code{nlar} methods.
}
\keyword{internal}
\keyword{ts}

