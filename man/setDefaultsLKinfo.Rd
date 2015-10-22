\name{setDefaultsLKinfo}
\alias{setDefaultsLKinfo}
\alias{setDefaultsLKinfo.default}
\alias{setDefaultsLKinfo.LKRectangle}
\alias{setDefaultsLKinfo.LKBox}
\alias{setDefaultsLKinfo.LKRing}
\alias{setDefaultsLKinfo.LKCylinder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Method for including default information in the LKinfo object.
%%  ~~function to do ... ~~
}
\description{This  method is used to define the various parts of the 
	LatticeKrig model by including specific components in the LKinfo object. 
	The specific instance of the function that is used depnds on the
	geometry that has been specified. Typically one would would not need
	to use this
    method for a new geometry. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
setDefaultsLKinfo(object, ...)
\method{setDefaultsLKinfo}{default}(object, ...)
\method{setDefaultsLKinfo}{LKRectangle}(object, ...)
\method{setDefaultsLKinfo}{LKBox}(object, ...)
\method{setDefaultsLKinfo}{LKRing}(object, ...)
\method{setDefaultsLKinfo}{LKCylinder}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  	An object of class LKinfo
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{Additional arguments to a specific method. 
%%     ~~Describe \code{\dots} here~~
}
}
\details{
	This method is used to include some default settings or components  in the
	 LKinfo  
	object and is called from \code{LKrigSetup}. 
	
	In creating a new geometry this method is not required but can be used to 
	include convenient default values and also some checks on the arguments passed to  LKrigSetup.
	
	For example with the LKRing geometry the defaults values follow the logic:
\itemize{	
	\item alpha set to  1.0 if not specified and nlevel is 1 
	
	\item a.wght set to 4.01 if not specified. 
	
     \item the fixed part of the model (fixedFunction) uses the function LKrigPeriodicFixedFunction which insures the first component is periodic.

     \item the transformation matrix, V has the default  diag( c(1,1)). Also in this
	function if V[1,1] is not 1 then an error is given.
	}
	 

}

\author{
	Doug Nychka
}

\seealso{
	\code{\link{LKrigSetup}}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial }
