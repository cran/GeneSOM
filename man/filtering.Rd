\name{filtering}
\alias{filtering}
\title{
Filter data before feeding som algorithm
}
\description{
Filtering data by certain floor, ceiling, max/min ratio, and max - min difference. 
}
\usage{
filtering(x, lt=20, ut=16000, mmr=3, mmd=200)
}
\arguments{
  \item{x}{
    a data frame or matrix of input data.
  }
  \item{lt}{
    floor value replaces those less than it with the value
  }
  \item{ut}{
    ceiling value replaced those greater than it with the value
  }
  \item{mmr}{
    the max/min ratio, rows with max/min < mmr will be removed
  }
  \item{mmd}{
    the max - min difference, rows with (max - min) < mmd will be removed
  }
}
\value{
  An dataframe or matrix after the filtering
}
\author{
  Jun Yan <jyan@stat.wisc.edu>
}
\seealso{
\code{\link{normalize}}.
}
\keyword{
  manip
}

