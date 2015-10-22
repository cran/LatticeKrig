\name{LKrigNormalizeBasis}
\alias{LKrigNormalizeBasis}
\alias{LKrigNormalizeBasisFast}                            
\alias{LKrigNormalizeBasisFast.default}
\alias{LKrigNormalizeBasisFast.LKRectangle}
\alias{LKRectangleSetupNormalization}

\title{Methods and functions to support normalizing to 
	marginal variance of one 

}
\description{The basis functions can be normalized so that the 
	marginal varaince of the process at each level and at all locations is
	one. A generic function \code{LKrigNormalizeBasis} computes this for
	any LatticeKrig model. However, in special cases the normalization can 
	be accelerated taking advantage of the geometry and the model for the
	lattice weights. This alternative is the method \code{LKrigNormalizeBasisFast}.

}
\usage{
LKrigNormalizeBasis( LKinfo, Level, PHI, ...)
LKrigNormalizeBasisFast(LKinfo, ...)
\method{LKrigNormalizeBasisFast}{default}(LKinfo, ...)
\method{LKrigNormalizeBasisFast}{LKRectangle}(LKinfo, Level, x,  ...)
LKRectangleSetupNormalization(mx, a.wght)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a.wght}{A.wght parameters.}
  \item{LKinfo}{ An LKinfo object.  NOTE:  Here  \code{choleskyMemory}, a Spam memory argument, can be a component of LKinfo and is subsequently passed
  	 through to the	(spam) sparse cholesky decomposition}
  \item{Level}{The multiresolution level.}
  \item{mx}{Matrix of lattice sizes.}
  \item{PHI}{Unnormalized basis functions evaluated at the locations 
  	to find the normalization weights.}
  \item{x}{Locations to find normalizaion weights.}
  \item{\dots}{Additional arguments for method. }
}
\details{
Normalization to unit variance is useful for reducing the artifacts of
the lattice points and creates a model that is closer to being stationary. 
The computation must be done for every point evaluated with the basis functions 
The basic calculation is 
\preformatted{
	tempB <- LKrigSAR(LKinfo, Level = Level)
    tempB <- LKrig.spind2spam(tempB)
    wght <- LKrig.quadraticform(t(tempB) %*% tempB, PHI = PHI, 
        choleskyMemory = choleskyMemory)
}
This generic method uses the low level utility \code{LKrig.quadraticform} that evaluates  diag( t(PHI) solve( Q.l) PHI ) where PHI are the basis functions evaluated at the locations and Q.l is the precision matrix for the lattice coefficients at level l. 
 

For constant a weights and for the rectangluar geometry one can use an eigendecomposition of the Kronecker product of the SAR matrix. 
}
\value{
A vector of weights. The basis are normalized by diving by the square root
of the weights.  
}
\author{
Doug Nychka
}

\keyword{ spatial }
