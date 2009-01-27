\name{MakeThSpec}
\alias{MakeThSpec}
\alias{makeThSpec}
\title{Specification of the threshold search}
\description{
This optional function allows the user to set different restrictions for the threshold grid search in function \code{\link{selectSETAR}}. 
}
\usage{
MakeThSpec(ngrid=c("All", "Half", "Third", "Quarter"), exact=NULL, int=c("from","to"), around="val",...)
}
\arguments{
\item{exact}{The user give an exact threshold value}
\item{int}{The user gives an interval to search inside}
\item{around}{The user gives an point to search around}
\item{ngrid}{The number of values to search for}
\item{...}{currently unused}
}
\value{
The input values are given as output after checking for consistency (only one of exact/int/around should be given). 
}
\details{
This function is just to check the inputs for the specification of the grid search. If not provided, the search will be in the biggest interval (ngrid ="All") between the minimum and maximum values. The user can reduce it by giving setting "Half" (only every two points is taken) and so on, or setting a number. 

The search can also be made around a point, or between two points. When between a point, the argument ngrid is still used, whereas for around, a value of 30 is taken as default value if ngrid is not specified by user. 
}
\seealso{\code{\link{selectSETAR}}}
\author{ Matthieu Stigler }
\examples{
#TODO
}
