\name{LKrig.coef}
\alias{LKrig.coef}
\alias{LKrig.lnPlike}
\alias{dfind2d}
\alias{dfind3d}
\title{
Internal functions for LKrig  and MLE.LKrig. 
}
\description{
Two internal functions for \code{LKrig} that estimate  the 
coefficients of the basis functions and compute the likelihood.
}
\usage{
LKrig.coef(Mc, wPHI, wT.matrix, wy, lambda)
LKrig.lnPlike(Mc, Q, y, lambda, residuals, weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{The ratio of the nugget variance (sigma squared) to the parameter
  controlling the marginal variance of the process (called rho in
  fields).} 

  \item{Mc}{Cholesky decomposition of regression matrix}
 
  \item{Q}{Precision matrix for coefficients.}

   \item{residuals}{Residuals from fitting spatial process.}

  \item{wPHI}{weighted matrix of basis functions. See \code{LKrig} for construction.}
  
  \item{wT.matrix}{weighted matrix of fixed part of estimate}

  \item{weights}{A vector that is proportional to the reciprocal
  variances of the errors. I.e. errors are assumed ot be uncorrealted
  with varainces sigma^2/weights }

  \item{wy}{weighted observations}

  \item{y}{Spatial observation}
}
\details{
\code{LKrig.coef} and \code{LKrig.lnPlike} are two low
level functions to find the basis function coefficients and the
evaluate the likelihood. \code{LKrig.lnPlike} also finds the
profile MLE of sigma and rho given a fixed value for lambda (and
alpha, a.wght).  See the code in LKrig and also
MLE.LKrig to see how these functions are used.
 
\code{dfind2d} is a fast FORTRAN suboutine to find nearest neighbors within a fixed distance and is called by \code{Wendland.basis}. \code{dfind3d} is currently not used but is intended for chordal distnace on the sphere. 

}
\value{
\code{LKrig.coef} a list with components d.coef the coefficients for covariates (Z)
and c.coef the basis function coefficients.
See \code{predict.LKrig} for how predicted values are found.

\code{LKrig.lnPlike} has the components:

 \code{lnProfileLike},  the log likelihood profiled for lambda, alpha and a.wght

 \code{rho.MLE}, the MLE of rho given   lambda, alpha and a.wght

 \code{shat.MLE}, the MLE of sigma  given   lambda, alpha and a.wght

 \code{quad.form}, the quadratic form in the exponent of the multivariate normal likelihood 

 \code{lnDetCov}, the log  determinant of the covariance matrix in the likelihood.

}
\author{
%%  ~~who you are~~
Doug Nychka
}
\seealso{
LKrig, LKrig.basis
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
