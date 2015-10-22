\name{LKinfoCheck}
\alias{LKinfoCheck}
\alias{LKinfoCheck.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Check the LKinfo object
%%  ~~function to do ... ~~
}
\description{
This method performs some simple checks on the LKinfo object before it is 
returned by LKrigSetup. One benefit is that it checks for the essential
components in each of the parts of the LKinfo object. 
Currently there is only a default method supplied. 

Here is an exmaple of how the checks work
\preformatted{
	LKinfo<- LKrigSetup( cbind(c(0,1)), LKGeometry="LKInterval",
    a.wght=5, nlevel=2, nu=1.5, NC=4)
    LKinfo$alpha
# corrupt it  
    LKinfo$alpha <- list( 1,.5,.25)
try( LKinfoCheck( LKinfo) )
}

}
\usage{
LKinfoCheck(object, ...)
\method{LKinfoCheck}{default}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  	An LKinfo object to be checked.
}
  \item{\dots}{
  	Additional arguments to method. 
}
}

\value{
There is no returned value. The benefit is the side effect of
an error message  and a "stop if the are problems in the LKinfo object. 
}

\author{
	Doug Nychka
}

\seealso{
\code{\link{LKrigSetup}}
}
\keyword{ spatial }
