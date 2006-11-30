\name{nlarDialog}
\alias{nlarDialog}
\title{GUI to nlar}
\description{
GUI interface to builtin NLAR models
}
\usage{
nlarDialog(series)
}
\arguments{
\item{series}{time series}
}
\details{
  Displays a GUI to \code{\link{nlar}}.
  Still under development. Is likely to change in future.
  Using the GUI, not all model options are available to the user.

  The finally fitted model is put in an object named \code{nlarModel} in the user workspace.
}
\value{
None.
}
\examples{
\dontrun{
if(interactive())
  nlarDialog(log(lynx))}
}
\author{Antonio, Fabio Di Narzo}
\keyword{ ts }
