\name{qerror}
\alias{qerror}
\title{ quantization accuracy }
\description{
  get the average distortion measure
}
\usage{
qerror(obj, err.radius=1)
}
\arguments{
  \item{obj}{ a `som' object }
  \item{err.radius}{ radius used calculating qerror }
}
\value{
  An average of the following quantity (weighted distance measure) over all x in the sample,
  \deqn{\sum ||x - m_i|| h_{ci}}{sum |x - m\_i| h\_{ci}}
  where \eqn{h_{ci}}{h\_{ci}} is the neighbourhood kernel for the ith code.
}
\author{
  Jun Yan <jyan@stat.wisc.edu>
}
\examples{
foo <- som(matrix(rnorm(1000), 100), 2, 4)
qerror(foo, 3)
}
\keyword{arith}
