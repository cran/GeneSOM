\name{som}
\alias{som.init}
\alias{som}
\alias{som.train}
\alias{som.update}
\alias{som.project}
\title{
Function to train a Self-Organizing Map
}
\description{
Produces an object of class "som" which is a Self-Organizing Map fit of the data.
}
\usage{
som.init(data, xdim, ydim, init="linear")
som(data, xdim, ydim, init="linear", alpha=NULL, alphaType="inverse",
neigh="gaussian", topol="rect", radius=NULL, rlen=NULL, err.radius=1,
inv.alp.c=NULL)
som.train(data, code, xdim, ydim, alpha=NULL, alphaType="inverse",
neigh="gaussian", topol="rect", radius=NULL, rlen=NULL, err.radius=1, inv.alp.c=NULL)
som.update(obj, alpha = NULL, radius = NULL, rlen = NULL, err.radius =
1, inv.alp.c = NULL)
som.project(obj, newdat)
}
\arguments{
  \item{obj}{a `som' object.}
  \item{newdat}{a new dataset needs to be projected onto the map.}
  \item{code}{a matrix of initial code vector in the map.}
  \item{data}{
    a data frame or matrix of input data.
  }
  \item{xdim}{
    an integer specifying the x-dimension of the map.
  }
  \item{ydim}{
    an integer specifying the y-dimension of the map.
  }
  \item{init}{
    a character string specifying the initializing method. The following
    are permitted:
    \code{"sample"} uses a radom sample from the data;
    \code{"random"} uses random draws from N(0,1);
    \code{"linear"} uses the linear grids upon the first two principle
    components directin.
  }
  \item{alpha}{
    a vector of initial learning rate parameter for the two training
    phases. Decreases linearly to zero during training.
  }
  \item{alphaType}{
    a character string specifying learning rate funciton type. Possible
    choices are linear function (\code{"linear"}) and inverse-time type
    function (\code{"inverse"}).
  }
  \item{neigh}{
    a character string specifying the neighborhood function type. The
    following are permitted:

    \code{"bubble"}
    \code{"gaussian"}
  }
  \item{topol}{
    a character string specifying the topology type when measuring
    distance in the map. The following are permitted:

    \code{"hexa"}
    \code{"rect"}
  }
  \item{radius}{
    a vector of initial radius of the training area in som-algorithm for
    the two training phases. Decreases linearly to one during training.
  }
  \item{rlen}{
    a vector of running length (number of steps) in the two training phases.
  }
  \item{err.radius}{
    a numeric value specifying the radius when calculating average
    distortion measure.
  }
  \item{inv.alp.c}{
    the constant C in the inverse learning rate function:
    alpha0 * C / (C + t);
  }
}
\value{
  `som.init' initializes a map and returns the code matrix. `som' does
  the two-step som training in a batch fashion and return a `som'
  object. `som.train' takes data, code, and traing parameters and
  perform the requested som training. `som.update' takes a `som' object
  and further train it with updated paramters. `som.project' projects
  new data onto the map.
  
  An object of class \code{"som"} representing the fit, which is a list
  containing the following components:
  
  \item{data}{the dataset on which som was applied.}
  \item{init}{a character string indicating the initializing method.}
  \item{xdim}{an integer specifying the x-dimension of the map.}
  \item{ydim}{an integer specifying the y-dimension of the map.}
  \item{code}{a metrix with nrow = xdim*ydim, each row corresponding to
    a code vector of a cell in the map. The mapping from cell coordinate
    (x, y) to the row index in the code matrix is:
    rownumber = x + y * xdim
  }
  \item{visual}{a data frame of three columns, with the same number of
    rows as in data: x and y are the coordinate of the corresponding
    observation in the map, and qerror is the quantization error
    computed as the squared distance (depends topol) between the observation
    vector and its coding vector.
  }
  \item{alpha0}{a vector of initial learning rate parameter for the two
    training phases.}
  \item{alpha}{a character string specifying learning rate funciton
    type.}
  \item{neigh}{a character string specifying the neighborhood function
    type.}
  \item{topol}{a character string specifying the topology type when measuring
    distance in the map.}
  \item{radius0}{a vector of initial radius of the training area in
    som-algorithm for the two training phases.}
  \item{rlen}{a vector of running length in the two training phases.}
  \item{qerror}{a numeric value of average distortion measure.}
  \item{code.sum}{a dataframe summaries the number of observations in
    each map cell.}
}
\references{
  Kohonen, Hynninen, Kangas, and Laaksonen (1995), SOM-PAK, the
  Self-Organizing Map Program Package (version 3.1).
  http://www.cis.hut.fi/research/papers/som\_tr96.ps.Z}
\author{
  Jun Yan <jyan@stat.wisc.edu>
}
\examples{
data(yeast)
yeast <- yeast[, -c(1, 11)]
yeast.f <- filtering(yeast)
yeast.f.n <- normalize(yeast.f)
foo <- som(yeast.f.n, xdim=5, ydim=6)
foo <- som(yeast.f.n, xdim=5, ydim=6, topol="hexa", neigh="gaussian")
}
plot(foo)
\keyword{cluster}
