\name{normalize}
\alias{normalize}
\title{
normalize data before feeding som algorithm
}
\description{
Normalize the data so that each row has mean 0 and variance 1.
}
\usage{
normalize(x, byrow=TRUE)
}
\arguments{
  \item{x}{
    a data frame or matrix of input data.
  }
  \item{byrow}{
    whether normalizing by row or by column, default is byrow.
  }
}
\value{
  An dataframe or matrix after the normalizing.
}
\author{
  Jun Yan <jyan@stat.wisc.edu>
}
\seealso{
  \code{\link{filtering}}.
}
\keyword{
  array, manip
}


