\name{tsDyn-package}
\alias{tsDyn-package}
\alias{tsDyn}
\title{Getting started with the tsDyn package}
\description{
Getting started with the tsDyn package
}
\details{
This package provide some tools inspired by nonlinear dynamics for the analysis-modelling of 
observed time series.\cr
For loading the package, type:\cr
\code{library(tsDyn)}\cr
A good place to start learning the package usage, is the vignette. It contains a more detailed guide on package contents, and an applied case study. At the R prompt, write:\cr
\code{vignette("tsDyn")}\cr
For a full list of functions exported by the package, type:\cr
\code{ls("package:tsDyn")}\cr
\cr
There is also an experimental GUI for built-in NLAR models. Call it with:\cr
\code{nlarDialog(timeSeries)}\cr
where \code{timeSeries} is an available time series object.

Each exported function has a corresponding man page (some man pages are in common to more functions). Display it by typing\cr
\code{help(functionName)}
}
\author{ Antonio, Fabio Di Narzo }
\seealso{
\code{\link{availableModels}} for listing all currently available NLAR models

\code{\link{autopairs}},\code{\link{autotriples}},\code{\link{autotriples.rgl}} for graphical explorative functions

\code{\link{llar}}, \code{\link{delta}}, \code{\link{delta.lin}} for nonlinearity checking tools
}
\keyword{ ts }
