\name{Wendland.basis}
\alias{Wendland.basis}
\alias{WendlandFunction}
\title{
Two dimensional radial basis functions based on a Wendland function.
}
\description{
Two dimensional radial basis functions based on a Wendland function and
using sparse matrix format to reduce the storage.
}
\usage{
Wendland.basis(x1, center, delta, max.points = NULL, mean.neighbor = 50, just.distance = FALSE)
WendlandFunction(d)
}

\arguments{
  \item{x1}{
   A 2 column matrix of   2-d locations to evaluate the basis functions. Each row of \code{x1} is a 
location. 

}
  \item{center}{
   A two column matrix specifying the basis function centers.

}
\item{d}{A vector of distances}
  \item{delta}{
   A vector of scale parameters for the basis functions  

}
  \item{max.points}{
   Maximum number of nonzero entries expected for the returned matrix.}

  \item{mean.neighbor}{
   Average number of centers that are within delta of each x1 location. For centers on a regular grid this
is often easy to estimate.}

  \item{just.distance}{
For debugging, don't evaluate the Wendland just return the pairwise distance matrix.}
}
\details{
This function finds the pairwsie distances between the points x1 and  centers and evaluates a Wendland 
function at these distances scaled by delta. In most applications delta is constant but a variable delta could be useful for 
lon/lat regular grids. The Wendland function is for 2 dimensions and smoothness order 2. See \code{WendlandFunction} for the 
polynomial form. This code has very simliar function to the fields function \code{wendland.cov}.  

In psudeo R code  for delta a scalar \code{Wendland.basis} evaluates as 
\preformatted{
  BigD<- rdist( x1,centers)
  WendlandFunction(BigD/delta)
}
The actual code uses a FORTRAN subroutine to search over distances less than delta and also returns the 
matrix in sparse format. 
}
\value{
For \code{Wendland.basis} a matrix in sparse format with number of rows equal to nrow(x1) and
columns equal to nrow(center).
}
\author{
Doug Nychka
}
\seealso{
LKrig.basis
}
\examples{
x<- cbind( runif(100), runif(100))
center<- expand.grid( seq( 0,1,,5), seq(0,1,,5))
# coerce to matrix
center<- as.matrix(center)

PHI<- Wendland.basis(x, center, delta=.5)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

