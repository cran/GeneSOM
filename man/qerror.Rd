\name{qerror}
\alias{qerror}
\title{ quantization accuracy }
\description{
  get the quantization accurarcy
}
\usage{
qerror(obj, radius=1)
}
\arguments{
  \item{obj}{ a `som' object }
  \item{radius}{ radius in calculating qerror }
}
\value{
  averaged quantizaiton error
}
\author{
  Jun Yan <jyan@stat.wisc.edu>
}
\examples{
foo <- som(matrix(rnorm(1000), 100), 2, 4)
qerror(foo, 3)
\keyword{arith}
