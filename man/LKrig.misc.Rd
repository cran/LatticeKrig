% # LatticeKrig  is a package for analysis of spatial data written for
% # the R software environment .
% # Copyright (C) 2012
% # University Corporation for Atmospheric Research (UCAR)
% # Contact: Douglas Nychka, nychka@ucar.edu,
% # National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
% #
% # This program is free software; you can redistribute it and/or modify
% # it under the terms of the GNU General Public License as published by
% # the Free Software Foundation; either version 2 of the License, or
% # (at your option) any later version.
% # This program is distributed in the hope that it will be useful,
% # but WITHOUT ANY WARRANTY; without even the implied warranty of
% # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% # GNU General Public License for more details.

% # Modified DH Dec 2013
% Earlier modifications DH 2/03/2012

\name{LKrig Miscellaneous Matrix Functions}
\alias{LKrig.rowshift.periodic}
\alias{ LKrig.shift.matrix}
\alias{LKrig.rowshift}
\alias{which.max.matrix}
\alias{which.max.image}
\alias{expandMList}
\alias{expandMatrix0}
\alias{expandMatrix}
\alias{repMatrix}
\title{
Miscellaneous internal functions for LatticeKrig package. 
}
\description{
Some utility functions used internally by higher level \code{LKrig}
functions. Currently these are simple functions that perform shifts of
a matrix.
}
\usage{
LKrig.shift.matrix( A, shift.row=0, shift.col=0, periodic=c(FALSE, FALSE))
LKrig.rowshift.periodic( A, shift.row)
LKrig.rowshift( A,shift.row,shift.col)
which.max.matrix(z)
which.max.image(obj)
expandMatrix0( A, B)
expandMatrix( ...)
expandMList( Mlist, byrow=TRUE)
}

\arguments{

\item{A}{A matrix.}

\item{byrow}{If TRUE matrices will be repeated row by row. If FALSE
this will be done column by column.}

\item{B}{Another matrix.}

\item{Mlist}{A list where each component is a matrix.}

\item{obj}{An image list with the usual components x, y, and z.}

\item{periodic}{ A vector of two logicals pertaining to rows and
columns. TRUE indicates an index where the shift will be periodic --
entries shifted beyond the dimensions will be wrapped to the other
side.}

\item{shift.row}{An integer that specifies the number of positions
 that the rows of the matrix are shifted.}

\item{shift.col}{An integer that specifies the number of positions
 that the columns of the matrix are shifted.}

\item{z}{A matrix.}

\item{\dots}{Matrices to be expanded.}
} \details{

\strong{Shift related:}
These functions are used to create the nearest neighbor indices for
the precision matrices.

\strong{Expand related:}
These functions are useful for creating a sets of covariance parametes that 
follow a factorial pattern. For example repeating the rows of the "alpha" paramters
as the "a.wght" parameters are varied. \code{expandMList} is particularly useful for creating
a factorial design of parameters to pass to LKrig.MLE for searching the likelihood.

}
\value{

\strong{Shift:}
A matrix of shifted values. Entries that are not defined due to
the shift are set to NA.  A column shift is done by a combination of
transpose operations and a row shift.

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

\strong{Expand:}
\code{ExpandMList}{Returns a list of matrices where the original matrices are repeated so that are
combinations of rows are represented. The example below illustrates.  \code{byrow=FALSE} does the repetition
by columns instead of rows.
}
\preformatted{
> A
     [,1] [,2]
[1,]    1    3
[2,]    2    4
> B
     [,1]
[1,]   11
[2,]   12
[3,]   13
> C
[1] 100
> expandMList( list( A=A, B=B, C=C))
$A
     [,1] [,2]
[1,]    1    3
[2,]    2    4
[3,]    1    3
[4,]    2    4
[5,]    1    3
[6,]    2    4

$B
     [,1]
[1,]   11
[2,]   11
[3,]   12
[4,]   12
[5,]   13
[6,]   13

$C
     [,1]
[1,]  100
[2,]  100
[3,]  100
[4,]  100
[5,]  100
[6,]  100
}
}

\author{
Doug Nychka
}


\keyword{spatial}
