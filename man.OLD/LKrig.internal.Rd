% corrections DMH 2/03/2012
\name{LKrig.coef}
\alias{LKrig.coef}
\alias{LKrig.lnPlike}
\alias{LKrig.traceA}
\alias{dfind2d}
\alias{dfind3d}
\title{
Internal functions for LatticeKrig package. 
}
\description{
Two internal functions for \code{LKrig} that estimate  the 
coefficients of the basis functions and compute the likelihood.
}
\usage{
LKrig.coef(Mc, wPHI, wT.matrix, wy, lambda, weights)
LKrig.lnPlike(Mc, Q, y, lambda, residuals, weights,
 sigma = NA, rho = NA)
LKrig.traceA( Mc, wPHI, wT.matrix, lambda, weights,NtrA, iseed=NA)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iseed}{Random seed used to generate the Monte Carlo samples. Keep the same to 
   compare results with mKrig and also for multiple values of lambda.}
  \item{lambda}{The ratio of the nugget variance (sigma squared) to the parameter
  controlling the marginal variance of the process (called rho in
  fields).} 

  \item{Mc}{Cholesky decomposition of regression matrix.}
  \item{NtrA}{Number of Monte  Carlo samples to estimate trace. Default is 20 in LKrig.}
  \item{Q}{Precision matrix for coefficients.}

   \item{residuals}{Residuals from fitting spatial process.}
  \item{rho}{If lambda is not specified the value of marginal varaince of the process.}
  \item{sigma}{If lambda is not specified the values of the measurement error standard deviation.}
  \item{wPHI}{Weighted matrix of basis functions. See \code{LKrig} for construction.}
  
  \item{wT.matrix}{Weighted matrix of fixed part of estimate.}

  \item{weights}{A vector that is proportional to the reciprocal
  variances of the errors. I.e. errors are assumed ot be uncorrelated
  with variances sigma^2/weights }

  \item{wy}{Weighted observations}

  \item{y}{Spatial observation}
}
\details{
\code{LKrig.coef} and \code{LKrig.lnPlike} are two low
level functions to find the basis function coefficients and to
evaluate the likelihood. The coefficents ( \code{c.mKrig}) are also found
because the provide short cut formulas for the standard errors and MLE esimates. 
These coefficients are identical to the basis coefficients (\code{c.coef}) found for usual Kriging in the mKrig function. \code{LKrig.lnPlike} also finds the
profile MLE of sigma and rho given a fixed value for lambda (and
alpha, a.wght).  See the code in LKrig and also
MLE.LKrig to see how these functions are used.

\code{LKrig.traceA} Finds an estimate of the effective degrees of freedom of the smoothing matrix based a simple Monte Carlo scheme. The smoothing matrix A is the matrix for fixed covariance 
parameters so that  y.hat = A y, where y.hat are the predicted values at the data locations. 
trace(A) is the effective degrees of freedom.  If e are iid N(0,1) then the expected value of 
t(e)\% * \% A \% * \% e is equal to the trace of A. This is the basis for estimating the  trace and the standard error for this estimate is based on \code{NtrA} independent samples. 
 
\code{dfind2d} is a fast FORTRAN subroutine to find nearest neighbors within a fixed distance and is called by \code{Wendland.basis}. \code{dfind3d} is currently not used but is intended for chordal distance on the sphere. 

}
\value{
\code{LKrig.coef} a list with components d.coef the coefficients for covariates (Z)
and c.coef the basis function coefficients.
See \code{predict.LKrig} for how predicted values are found.

\code{LKrig.lnPlike} has the components:

 \code{lnProfileLike},  the log likelihood profiled for lambda, alpha and a.wght

 \code{rho.MLE}, the MLE of rho given lambda, alpha and a.wght

 \code{shat.MLE}, the MLE of sigma given lambda, alpha and a.wght

 \code{quad.form}, the quadratic form in the exponent of the multivariate normal likelihood 

 \code{lnDetCov}, the log determinant of the covariance matrix in the likelihood.

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
