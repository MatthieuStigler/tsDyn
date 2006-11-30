\name{toLatex.setar}
\alias{toLatex.setar}
\title{Latex representation of fitted setar models}
\description{
Latex representation of fitted setar models
}
\usage{
\method{toLatex}{setar}(object, digits=3, \dots)
}
\arguments{
  \item{object}{fitted setar model (using \code{\link{nlar}})}
  \item{digits, ...}{options to be passed to \code{\link{format}} for
    formatting numbers}
}
\seealso{
\code{\link{setar}}, \code{\link{nlar-methods}}
}
\author{ Antonio, Fabio Di Narzo }
\examples{
mod.setar <- setar(log10(lynx), m=2, thDelay=1, th=3.25)
toLatex(mod.setar)
}
\keyword{ ts }
