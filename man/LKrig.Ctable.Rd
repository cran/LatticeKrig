\name{LKrig.Ctable.10.3}
\alias{LKrig.Ctable.10.3}
\alias{LKrig.Ctable.10.4}
\docType{data}
\title{Table of LKrig parameters to approximate the Matern covariance }
\description{
Tables specifying the paramters alpha, a.wght and 
}
\usage{
data(LKrig.Ctable.10.3)
data(LKrig.Ctable.10.4)
 }
\format{
  The table format is a list with components:
"NC"         "NL"         "range"      "smoothness" "kappa"     
"a.wght"     "alpha"      "par.grid"   "par.glist" 

}
\details{
The parameters are found by least squares fitting to the Matern covariance for a grid of ranges.
These tables are used by the function \code{Matern2LKrig} to determine an approximate 
LKrig covariance. The naming convention for these tables is LKrig.Ctable.X.Y with X
being the number of basis functions in first level and Y being the number of levels.
Note that both a.wght and kappa are equivalent parameters. 

Currently the approximation is only implemented for the exponential (smoothness=.5).
 
}


\keyword{datasets}
