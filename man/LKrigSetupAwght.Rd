\name{LKrigSetupAwght}
\alias{LKrigSetupAwght}
\alias{LKrigSetupAwght.default}
\alias{LKrigSetupAwght.LKRectangle}

\title{Method to create a.wght component from the \code{LKinfo}
	object.
%%  ~~function to do ... ~~
}
\description{This method takes a vector or more complicated object and
	based on the geometry creates a list with the a.wght information.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LKrigSetupAwght(object, ...)
\method{LKrigSetupAwght}{default}(object, ...)
\method{LKrigSetupAwght}{LKRectangle}(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  		The partial LKinfo object created within \code{LKrigSetup}
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{
  	Any additional arguments to this method
%%     ~~Describe \code{\dots} here~~
}
}
\details{
	The simplest function of this method is to convert the sequence of 
 	a.wght values into a list. If only a scalar a.wght value is supplied then the default method just repeats this for each level. The attribute
 	\code{fastNormalize} (either TRUE or FALSE) is attached to this 
 	list to indicate how the marginal variance of the process should be found.
 	  
\preformatted{ 	 	
  LKinfo<- LKrigSetup( x,LKGeometry="LKInterval", alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 3)
  LKrigSetupAwght( LKinfo)
  
[[1]]
[1] 4.5

[[2]]
[1] 4.5

[[3]]
[1] 4.5

attr(,"fastNormalize")
[1] FALSE
}                   
 	Currently the only geometry with fastNormalization being TRUE is for a
 	rectangular domain. 
 	
 	 For the LKRectangle geometry, however, more complicated
	anisotropic and nonstationary a.wght specifications are possible. 
	See \code{\link{LKrig}} for details. Also in the case that the
    fastNormalization is TRUE for rectangles several more attributes are added to the  \code{a.wght} list 	that precompute some matrices of the SAR. 
	
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A list with \code{nlevel} components. The attribute \code{fastNormalize}
	is added to this list. In the case that the geometry is LKRectangle several
	more attributes are added indicating the type of covariance model and
	possibly an eigen decomposition of the SAR matrix exploiting Kronecker
	products.  
}

\author{
	Doug Nychka
%%  ~~who you are~~
}


\seealso{
	\code{\link{LKrigSetup}}, 	\code{\link{LKrigSetupAlpha}},
	\code{\link{LKrigSAR}},  \code{\link{LKrig}}

}
\examples{
  x<- cbind( c(0,1))
  LKinfo<- LKrigSetup( x,LKGeometry="LKInterval", alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 3)
  a.wghtList<- LKrigSetupAwght( LKinfo)
  
  x<- cbind( c(0,1), c(0,1))
  LKinfo<- LKrigSetup( x, alpha=c( 1,.2,.01),
                   nlevel=3, a.wght=4.5, NC= 3)
  a.wghtList<- LKrigSetupAwght( LKinfo)
# see   
  names(attributes( a.wghtList))
 
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }

