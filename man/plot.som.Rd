\name{plot.som}
\alias{plot.som}
\title{ Visualizing a SOM }
\description{
 Plot the SOM in a 2-dim map with means and sd bars.
}
\usage{
\method{plot}{som}(obj, sdbar=1, ylim=c(-5, 5), color=TRUE,
ntik=3, yadj=0.1, xlab="", ylab="", ...)
}
\arguments{
  \item{obj}{ a som object }
  \item{sdbar}{ the length of sdbar in sd, no sdbar if sdbar=0 }
  \item{ylim}{ the range of y axies in each cell of the map }
  \item{color}{ whether or not use color plotting }
  \item{ntik}{ the number of tiks of the vertical axis }
  \item{yadj}{ the proportion used to put the number of obs }
  \item{xlab}{ x label }
  \item{ylab}{ y label }
  \item{\dots}{ other options to plot }
}
\author{ Jun Yan <jyan@stat.wisc.edu> }
\examples{
foo <- som(matrix(rnorm(1000), 250), 3, 5)
plot(foo, ylim=c(-1, 1))
}
\keyword{
  hplot
}