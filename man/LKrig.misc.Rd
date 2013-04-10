% corrections DMH 2/03/2012
\name{LKrig Miscellaneous Functions}
\alias{LKrig.rowshift.periodic}
\alias{ LKrig.shift.matrix}
\alias{LKrig.rowshift}
\title{
Internal functions for LatticeKrig package. 
}
\description{
Some utility functions used internally by 
higher level \code{LKrig} functions. Currently these are simple functions
that perform shifts of a matrix.
}
\usage{
  LKrig.shift.matrix( A, shift.row=0, shift.col=0, periodic=c(FALSE, FALSE))
  LKrig.rowshift.periodic( A, shift.row)
  LKrig.rowshift( A,shift.row,shift.col)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{A}{A matrix.}
\item{periodic}{ A vector of two logicals pertaining to rows and columns. TRUE indicates an index
where the shift will be periodic -- entries shifted beyond the dimensions will be wrapped to the
other side.} 
\item{shift.row}{An integer that specifies the number of positions that the rows of the matrix
 are shifted.}
\item{shift.col}{An integer that specifies the number of positions that the columns of the matrix
 are shifted.}
}
\details{
These functions are used to create the nearest neighbor indices for the precision matrices. 
}
\value{
A matrix of shifted values. For entries that are not defined due to the shift are set to NA. 
A column shift is done by a combination of transpose operations and a row shift. 
\preformatted{
A<- matrix( 1:12,3,4)
A
     [,1] [,2] [,3] [,4]
[1,]    1    4    7   10
[2,]    2    5    8   11
[3,]    3    6    9   12

#shift of 2 for rows:
 LKrig.rowshift( A, 2)
    [,1] [,2] [,3] [,4]
[1,]   NA   NA   NA   NA
[2,]   NA   NA   NA   NA
[3,]    1    4    7   10

#periodic case
LKrig.rowshift.periodic( A, 2)
     [,1] [,2] [,3] [,4]
[1,]    2    5    8   11
[2,]    3    6    9   12
[3,]    1    4    7   10
}

}
\author{
%%  ~~who you are~~
Doug Nychka
}
\seealso{
LKrig.MRF.precision
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
