% corrections DMH 2/03/2012
\name{LatticeKrig}
\alias{LatticeKrig}
\alias{LKrig}
\alias{Matern2LKrig}
\alias{predict.LKrig}
\alias{predict.se.LKrig}
\alias{print.LKrig}
\alias{print.LKinfo}
\alias{surface.LKrig}
\alias{LKrig.predict.surface}
\title{Spatial prediction and inference using a compactly supported
multi-resolution basis and a lattice model for the basis coefficients.
}
\description{
A variation of Kriging with fixed basis functions that uses the Wendland covariance kernel to 
create a regular set of basis functions on a grid. The coefficients of these
basis functions are modeled as a Markov random field (MRF).  
}
\usage{
LKrig( x,y=NULL, weights = rep(1, nrow(x)),Z=NULL,
                        LKinfo=NULL, 
                        iseed=123,NtrA=20,
                        use.cholesky=NULL, return.cholesky=TRUE,
                        NC, nlevel, a.wght, alpha, lambda=NA, sigma=NA,
                        rho=NA,  rho.object=NULL,
                        overlap=2.5, normalize=TRUE,edge=TRUE,
                        RadialBasisFunction="WendlandFunction",
                        verbose=FALSE)

Matern2LKrig( range, smoothness,x = NULL, xlim=c(-1,1), ylim=c(-1,1), nlevel, NC=10,
                                   check=FALSE, NP=500, NCtable=10)

\method{predict}{LKrig}( object, xnew=NULL,Znew = NULL, drop.Z=FALSE, ...)

\method{predict.se}{LKrig}( object, xnew=NULL, Znew = NULL, drop.Z = FALSE, verbose = FALSE, ...)

\method{surface}{LKrig}( object, ...)

LKrig.predict.surface( object, ...)

\method{print}{LKrig}( x, digits=4, ...)
}
\arguments{
  
  \item{alpha}{Weights of precision matrices for each level of resolution. This can be a vector or a list with \code{nlevel} components specifying how the alpha varying spatially. See details below. }

  \item{a.wght}{The "a" spatial autoregressive parameter for a first order Markov Random field
 (see details  below) . This controls the spatial dependence and must be greater than or equal to 4. For \code{a.wght} = 4 normalize should be FALSE.
  If there are multiple levels this can be a vector. See Details below.}



  \item{digits}{Number of digits in printed output.}
   
  \item{drop.Z}{If true the fixed part will only be evaluated at the
spatial drift polynomial part of the fixed model. The contribution from the other
covariates ( e.g elevation) in the fixed component will be omitted.}

  \item{edge}{If TRUE an adjustment is made in precision matrices for edge effects.}

%  \item{grid.info}{A list with information to construct the regular grid for the basis function centers.
%  The components are: \code{xmin, xmax, ymin,ymax, delta}. If not passed this list is created by
%  \code{LKrig.setup} from \code{x} and \code{NC}.}

 \item{iseed}{Random seed used in the Monte Carlo technique
  for approximating the effective degrees of freedom (trace of the
  smoothing matrix) for the GCV criterion. If NA, no seed is set. }

 \item{lambda}{The ratio of the nugget variance (called sigma squared in \code{fields} and \code{LatticeKrig}) to the parameter 
  controlling the marginal variance of the process (called rho in
  \code{fields} and \code{LatticeKrig}). If sigma and rho are both specified then \code{lambda == sigma^2/rho}.} 

 \item{LKinfo}{A list that is an alternate form of specifying the multiresolution basis. This is usually 
        created by the function \code{LKrig.setup}. If NULL, this list is created and returned as 
a component of the LKrig object. }

   \item{NC}{Maximum number of grid points in one dimension. If the
domain is square there will be NC*NC basis functions/ lattice points in the
coarsest level of resolution. If domain is rectangular, the smaller dimension will have less than NC
points.}

  \item{nlevel}{Number of levels for the multiresolution basis. Each level increases the
number of basis functions by a factor of 4.}

  \item{normalize}{If \code{TRUE} basis functions are normalized so that the 
    marginal variance of the process covariance is constant and equal to rho. This
    normalization avoids some edge and wavy artifacts from using a discrete
    set of basis functions.}

  \item{NtrA}{Number of random samples used in Monte Carlo method for
determining effective degrees of freedom. }
  
  \item{object}{The \code{LKrig} object.}

  \item{overlap}{The overlap between basis functions. This scaling is based on centers being spaced 
1 unit distance apart and the Wendland function decreasing to zero at 1 unit distance. A scaling/overlap of 2.5 (the default)
implies that the support of the basis functions will extend to a disc of radius 2.5.
 We recommend that this parameter not be changed unless one is sure of the effect on the 
implied spatial covariance.}

  \item{return.cholesky}{If TRUE the cholesky decomposition is included in the output list
    (with the name \code{Mc}).  Set to FALSE to avoid much larger objects when the decomposition is not needed f     or later computations (e.g. MLE search or standard errors). This option is often paired with a
     subsequent call to
     \code{LKrig} with \code{use.cholesky}. See the MLE.LKrig code for details. }
  \item{RadialBasisFunction}{An R function that defines the radial basis function to generate the multiresolution basis. Is assumed that this
            function is compacted supported on the interval [0,1]. See WendlandFunction for an example.}
  \item{rho}{Value of rho (sill variance) to use for spatial estimate. If omitted this is estimated as the MLE based on the value of lambda.}
  \item{rho.object}{A prediction object that defines a spatially varying component for the 
           marginal variance of the process. The object should be such that the code 
   \code{predict(rho.object, x1)} will evaluate the function at 2-d locations \code{x1}. See details  below.
} 
  \item{sigma}{Value of sigma (nugget variance) to use for spatial estimate.
            If omitted this is estimated as the MLE based on the value of lambda.}
  \item{verbose}{If \code{TRUE} print out results as they are computed in loop.}
  
  \item{use.cholesky}{Use the symbolic part of the Cholesky decomposition passed in  this argument.}

  \item{weights}{A vector that is proportional to the reciprocal
  variances of the errors. I.e. errors are assumed to be uncorrelated
  with variances sigma^2/weights. }

  \item{xnew}{Matrix of locations for prediction.}
 
  \item{x}{Spatial locations of observations.}
 
  \item{y}{Spatial observations.}
 
  \item{Z}{Linear covariates to be included in fixed part of the model
  that are distinct from the default first order polynomial in \code{x} (i.e. the spatial drift).}
  \item{Znew}{Values of covariates, distinct from the spatial drift for predictions off of data locations}

  \item{\dots}{Additional arguments to pass to generic methods. See help(predict.surface)
    for the details in calling LKrig.predict.surface}

  \item{range}{Range parameter for Matern target covariance in the LKrig approximation}
  \item{smoothness}{Matern smoothness parameter. Currently is only .5 (exponential).}
  \item{xlim}{X limits of domain}
  \item{ylim}{Y limits of domain}
  \item{check}{If TRUE returns the LKrig approximate covariance and the target 
              Matern on a grid of distances from the center of domain}
  \item{NP}{Number of points to evaluate the LKrig and Matern covariances}
  \item{NCtable}{Value of NC that is used to create table of approximations. Currently 
                  only NCtable=10 is supported.}
}
\details{
This method combines compactly supported basis
functions and a Markov random field covariance model to provide
spatial analysis for very large data sets. The key is to make sure all components are sparse matrices where linear algebra is required for the computations. The model for the spatial
field (or spatial process) is
 
 f(x) = P(x) + Z d + sum Phi.j(x) c.j

and the data model is

 Y.k = f(x.k) + e.k

with e.k uncorrelated normal errors with variance sigma^2/weights.
Where P is a first order linear polynomial, Z a matrix of covariates,
d the fixed coefficients, Phi.j are fixed basis functions and c.j are
random coefficients.

\strong{Radial basis functions:}
 The basis
 functions are two dimensional radial basis functions (RBF) proportional to 
 the Wendland covariance. There is one important difference in the basis construction
 -- a normalization -- and this aspect makes it different from a simple radial basis function specification and this is  described below. 
 For a single level (\code{nlevel==1}) the
 RBFs are centered at the grid points and with radial support
 \code{delta*overlap} where \code{delta} is the spacing between grid
 points and overlap has the default value of 2.5. For the multi-resolution
 each subsequent level is based on a grid with delta divided by
 2 (see example below) .

\strong{Markov random field (MRF) for the coefficients:}
The prior on c is a multivariate normal, with a mean of zero and the
inverse of the covariance matrix, Q, also known as the precision
matrix, is assumed to be block diagonal corresponding to the
organization of the basis functions into levels of resolution.  If
\code{nlevels} are specified, the ith block is a precision matrix
based on a first order spatial autoregression with \code{a.wght[i]}
being related to the autoregressive parameter. The specific precision matrix for
each block (level), Q.i, is implemented by
\code{LK.MRF.precision}. Briefly, this matrix has the form Q.i=
t(B)\%*\%B and is multiplied by the relative weight \code{1/alpha[i]}.
In the case of a spatially varying \code{alpha} the modification is

Q.i = t( diag( 1/sqrt(alpha[[i]])\%*\%B) \%*\% diag( 1/sqrt(alpha[[i]])\%*\%B

Note that the spatial variation of alpha is parameterized for the node locations of the 
basis functions, i.e. the locations of the lattice points for the GMRF. Of course if alpha
is constant for each  level this more complicated formula for Q.i reduces to 
 \code{1/alpha[i] * t(B)\%*\%B }. 

If the element in  ith precision/covariance matrix  has a marginal variance of say V
 then \code{alpha[i]}V will be the marginal variance of the ith component. 
More will be written about parametrizing the marginal variances below. 
Each row of B, corresponding to a point in the lattice in the
interior, is "a" (\code{a.wght[i]}) on the diagonal and  -1 at each of the
four nearest neighbors of the lattice point. Points on the edges just
have less neighbors if \code{edge} is set to FALSE.
Following the format in Lindgren and Rue (2010)  let a= 4 + E with E greater than or equal to  0.
(Where code{sqrt(1/E)} can be associated with an approximate range parameter.)
Some schematics for filling in the B matrix are:

\tabular{cll}{

       -1        \tab         -1           \tab             (4+E)  -1 \cr
                 \tab                      \tab                      \cr
  -1   (4+E)  -1  \tab        (4+E)  -1      \tab            -1        \cr
                 \tab                      \tab                      \cr
       -1        \tab        -1            \tab                      \cr
                 \tab                      \tab                      \cr
Interior point   \tab     Left edge when   \tab         Upper left corner when  \cr
                 \tab     edge is FALSE       \tab         edge is FALSE
}
 

If \code{edge} is TRUE the edges are weighted to reflect other boundary conditions. 

\tabular{cll}{

       -1        \tab         -.5           \tab             (1+E/4)  -.5 \cr
                 \tab                      \tab                      \cr
  -1   (4+E)  -1  \tab        (2+E)  -1      \tab            -.5        \cr
                 \tab                      \tab                      \cr
       -1        \tab        -.5            \tab                      \cr
                 \tab                      \tab                      \cr
Interior point   \tab     Left edge when   \tab         Upper left corner when  \cr
                 \tab     edge is TRUE       \tab         edge is TRUE
}
 
 
 The assumption is that for the coefficients associated with the ith
block of basis functions (c.i), t(B)\%*\%c.i will be uncorrelated.
This description is a spatial autoregressive model (SAR). The matrix Q
will of course have more nonzero values than B and the entries in Q
can be identified as the weights for a conditional autoregressive
model (CAR). Moreover, the CAR specification defines the neighborhood
such that the Markov property holds. Another way of viewing this model
is in terms of the precision and it is an elementary property of the
multivariate normal that relates the nonzero elements of Q to
conditional independence of components of the multivariate random
vector.  Values for \code{a.wght[i]} that are greater than 4 give
reasonable covariance models.  Moreover setting a.wght to 4 and
normalizing to FALSE in the call to LKrig will give a thin-plate
spline type model that does not have a range parameter. An approximate
strategy is to set a.wght close to 4 and get some benefit from the
normalization for edge effects.


Finally it should be noted that the specification of the basis and MRF
is through the components of the list LKinfo, a required component for
many LatticeKrig functions. High level functions such as \code{LKrig}
still require a minimal amount of information and combined with
default choices create the LKinfo list. One direct way to create the
complete list is to use \code{LKrig.setup} as in the example below.
For \code{nlevel==1} one needs to specify \code{a.wght}, \code{NC},
and also \code{lambda} related to the measurement error variance.  For
a multiresolution setup, besides these parameters, one should consider
different values of \code{alpha} keeping in mind that if this vector
is not constrained in some way ( e.g. \code{sum(alpha)==1}) it will
not be identifiable from \code{lambda}.

The overall design of this function does not require that the centers
actually define a grid, but then the construction of the matrix B may
need to be modified in the internal function \code{LK.MRF.precision}
and the basis function centers need to be include in the \code{LKinfo}
list.

 \strong{The covariance derived from the MRF and basis functions:} An
important part of this method is that the MRF describes the
coefficients of the basis functions rather than the field itself. Thus
in order to determine the covariance for the observed data one needs
to take in account both the MRF covariance and the expansion using the
basis functions. The reader is referred to the function \code{LKrig.cov}
for an explicit code that computes the implied covariance function for
the process f.  Formally if P is the covariance matrix (the inverse of
Q) for the coefficients then the covariance between the field at two
locations x1 and x2, will be

   sum_ij P_ij Phi.i(x1) Phi.j(x2)

 Moreover under the assumption that coefficients at different levels
are independent this sum can be decomposed into sums over the separate
levels.  The function \code{LKrig.cov} evaluates this formula based on
the LKrig object (\code{LKinfo}) at arbitrary groups of locations
returning a cross covariance matrix. \code{LKrig.cov.plot} is a handy
function for evaluating the covariance in the x and y directions to
examine its shape.  The function \code{LK.cov} is also in the form to
be used with conventional Kriging codes in the fields package (loaded
by LatticeKrig) \code{mKrig} or \code{Krig} and can be used for
checking and examining the implied covariance function.


The unnormalized basis functions result in a covariance that has some
non-stationary artifacts (see example below) For a covariance matrix P
and for any location x one can evaluate the marginal variance of the
unnormalized basis functions for each multiresolution level.  Based
this computation there is a weighting function, say w.i(x) so that
when the normalized basis w.i(x) Phi.i(x) is used the marginal
variance for the multiresolution process at each level will be one.
This makes the basis functions dependent on the choice of Q and
results in some extra overhead for computation.  But we believe it is
useful to avoid obvious artifacts resulting from the use of a finite
spatial domain (edge effects) and the discretization used in the basis
function model. This is an explicit way to make the model stationary
in the marginal variance with the result that the covariance also
tends to be closer to a stationary model. In this way the
discretization and edges effects associated with the MRF can be
significantly diminished.

The default in \code{LKrig} is \code{normalize = TRUE}. It is an open
question as to whether all levels of the multi-resolution need this
kind of explicit normalization. There is the opportunity within the
\code{LKrig.basis} function to only normalize specific levels with the
\code{normalize} being extended from a single logical to a vector of
logicals for the different levels.
To examine what these edge effect artifacts are like 
the marginal variances for a 6 by 6 basis is included at the end of the Examples Section.

Given that the process at each level has been normalized to have
marginal variance of one there are several other points where the
variance can be modified. The variance at each level is scaled by the
vector of parameters \code{alpha} and the marginal variance of the
process is scaled by \code{rho}. There is then a problem of
identifiability between these and a good choice is to constrain
\code{sum(alpha) ==1} so that \code{rho} is associated with the
marginal variance of the full spatial process. A final option in
specifying the marginal variance is to prescribe a spatially varying
multiplier. This component is specified by the object
\code{rho.object}. By default this is not included (or assumed to be
identically one) but, if used, the full specification for the marginal
variance of the spatial process at location \code{x} is formally:

\code{rho * predict(rho.object, x) * sum( alpha)}

provided the basis functions are normalized. 

To keep the documentation modular, help on evaluating the LKrig covariance, unconditional and conditional simulation can be found in help(LKrig.cov).

\strong{The smoothing parameter lambda and effective degrees of freedom}
Consistent with other fields package functions, the two main
parameters in the model, sigma^2 and rho are parameterized as lambda
= sigma^2/rho and rho. The MLEs for rho and sigma can be written in closed form as a function 
of lambda.  These estimates can be substituted into the full likelihood to give a profiled version that can numerically be maximized over lambda. See help( MLE.search.LKrig) for some examples of this
operation. lambda is best varied on a log scale but is unfortunately  difficult to interpret independent of the particular set of locations, sample size and covariance. A more useful interpretation of lambda is through the effective degrees of freedom and this value is found by default using a Monte Carlo technique. The effective degrees of freedom will vary with the dimension of the fixed regression part of the spatial model ( typical 3 = constant + linear terms) and the total number of observations. It can be interpreted as the approximate number of "parameters" needed to represent the spatial prediction surface. For a fixed covariance model the spatial estimate at the observation locations can be 
represented as  f hat = A(lambda) y where A is a matrix that does not depend on observations. 
The effective number of degrees of freedom is the trace of A(lambda) in analogy to the least squares regression "hat" matrix and has an inverse, monotonic relationship to lambda. The Monte Carlo method
uses the fact that if e are iid N(0,1)  E( t(e) A(lambda) e) = trace( A(lambda)). 


\strong{Spatial prediction:} 
Given this model, it is possible to compute the conditional expectation of f
given Y and also evaluate the profile likelihood for lambda, alpha, 
a.wght  and obj.rho (if appropriate). Because both the basis functions and
Q are sparse, it is possible to compute the estimates and likelihood for
large numbers of spatial locations. This model has an advantage over covariance tapering or
compactly supported covariance functions (e.g. fastTps), because the
implied covariance functions can have longer range correlations.

\strong{Descriptions of specific functions and objects:}

\strong{\code{LKrig}: } Find spatial process estimate for fixed covariance specificed by 
\code{ nlevel, alpha, a.wght, NC,  and lambda}  or this information in an \code{LKinfo} list. 

\strong{\code{LKinfo}:} This object is a "list of lists" and collects all the information to describe the covarince model and the observatonal model. In most case one would supply some of the important components of this andthe functions \code{LKrig.setup} and \code{LKrig} will fill in  the rest with their defaults. The basic design is to pass the LKinfo object to subsequent functions to avoid keeping track of many parts of the covariance specification. For example simulating an unconditional sample of the process is specified by
 passing the LKinfo object to 
\code{LKrig.cov}. A print method for an LKrig object prints a summary without listing large grids and 
parameter matrices. 
The current components of \code{LKinfo} are:

\describe{
\item{mx my}{ The grid for the basis functions is mx X my}
\item{nlevel}{Number of multiresolution levels}
\item{delta}{Spacing between grid point centers}
\item{m}{Total number of basis functions}
\item{offset}{Indices for coefficients by level. offset[k]+1 is the position of the first coefficient in level k}
\item{grid.info}{A list with min and max values for x and y grids and the spacing.}
\item{grid}{ A list with nlevel components:  grid[[k]] is a list with components  \code{x} and \code{y} specifying the grid for the basis function centers at level k.}
\item{overlap}{Amount the basis functions overlap.}
\item{alpha}{A list with nlevel components with alpha[[k]] being the scalar or matrix giving the values for alphat the kth level.}
\item{a.wght}{A list with nlevel components with a.wght[[k]] being the scalar or matrix giving the values for alphat the kth level.}
\item{lambda}{Value of lambda. A default values of  NA is OK if it has not been estiamted or set.}
\item{ sigma rho}{Values of the covariance parameters for measurement standard deviation and the process marginal variance. These can be NA if not yet determined.}
\item{normalize}{A logical value with TRUE being to normalize the basis function to have constant marginal variance}
\item{ normalize.level}{ A logical vector of length nlevel to determine which levels of basis functions should be normalized.}
\item{edge}{A logical if TRUE edge corrections are made to precision matrix.}
\item{scalar.alpha}{Logical value with TRUE indicating that alpha is a scalar at each level -- i.e. not a matrix.}
\item{scale.basis}{Logical if TRUE indicates rho is a prediction obj defining a surface.}
\item{rho.object}{The prediction object that define the value of rho at every point in the domain according to
the predict function.}
\item{RadialBasisFunction}{Name of radial basis function used as template for generating basis functions.}
}


\strong{\code{predict.LKrig, predict.se.LKrig}: } 
These functions evaluate the model at the 
the data locations or at \code{xnew} if it is included. Note the use of the  \code{drop.Z} argument to either include the covariates or just restrict the computation to the spatial 
drift and the smooth component. If \code{drop.Z} is FALSE then then \code{Znew} needs to be included for predictions off of the observation locations. The standard errors are computed under the assumption that the covariance is known, that it is the TRUE covariance for the process, and  that the process and measurement errors are multivariate normal. The formula that is used involves some recondite shortcuts for efficiency but has been checked against the standard errors found from an alternative formula in the fields Krig function. (See the script Lkrig.se.tests.R in tests for details.)


\strong{\code{Matern2LKrig}: } This is an experimental function to approximate a Matern covariance
function with a multiresolution, LKrig type covariance. Currently only tables for 
the exponential and 3 and 4 levels have been computed and \code{NC} must be 10 or larger. Also this is computed for a limited set of range parameters. An error is returned if the range is larger or smaller than the limits covered by the tables. 
Despite all these restrictions the function is useful. To check the quality of the approximation if \code{check=TRUE} the approximate LKrig covariance and the 
target Matern evaluated on a grid of distances are returned.
For example to test the approximation to a range of .2 on the domain [-1,1]X[-1,1]
(the default) using NC=10 and 4 levels. 
\preformatted{
  cov.approx<- Matern2LKrig( range=.2, smoothness=.5, NC=10, nlevel=4, check=TRUE)
  matplot( cov.approx$d, cov.approx$y, type="l")
}
Dashed line is LKrig approximation and solid is Matern with distance evaluated in the
X direction. See the components \code{d2} and \code{y2} to check the fit in the 
Y direction. Note that the implied covariance for the LKrig model can be evaluated at
arbitrary points using the \code{LKrig.cov} function. 


}
\value{
\strong{\code{LKrig: }}  
An LKrig class object with components for evaluating the estimate at arbitrary locations, 
describing the fit and as an option ( with \code{Mc.return=TRUE}) the Cholesky decomposition to allow for fast updating with new values of lambda, alpha and a.wght. 
The "symbolic" first step in the sparse Cholesky decomposition can also be used to compute the sparse Cholesky decomposition for a different positive definite 
matrix that has the same pattern of zeroes. This option is useful in computing the likelihood
under different covariance parameters. 
The returned component \code{LKinfo} is a list with the information that describes the layout of the multiresolution basis functions and the covariance parameters for the MRF. 
(see \code{LK.basis} as an example.)

\strong{\code{predict.LKrig, predict.se.LKrig: }} A vector of predictions or standard errors. 

\strong{\code{LKrig.predict.surface}: } A list in image format (components \code{x,y,z}) of the surface evaluated on a regular grid. This surface can then be plotted using 
several different core and fields functions e.g. \code{image}, \code{image.plot}\code{contour}, \code{persp}, \code{drape.plot}. The \code{surface} method just calls this
function and then a combination of the image and contour plotting functions. 



}
\author{
Doug Nychka
}

\seealso{
mKrig, Krig, fastTps, Wendland, LKrig.coef, Lkrig.lnPlike, LK.MRF.precision, LK.precision
}
\examples{
# Load ozone data set
  data(ozone2)  
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
# Find location that are not 'NA'.
# (LKrig is not set up to handle missing observations.)
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]

##############################################################################################
# Predict ozone concentrations over the domain using a single level lattice model
# If the domain is square there will be 30X30 = 900 basis functions and the
# a.wght is 5, smoothing parameter lambda = sigma^2 / rho is set to .01
##############################################################################################

  obj<- LKrig(x,y,NC=20,nlevel=1, alpha=1, lambda=.01, a.wght=5)
# summary of fit
  print(obj)
# Plot fitted surface.
  surface(obj) # see also predict.surface to just evaluate on a grid
  US(add=TRUE)
# Check effective degrees of freedom in surface.
  obj$trA.est
# NOTE: obj$LKinfo  has the full description of the Lattice Krig covariance model.
# A cleaner way to fit the same model separating creating the LKinfo object and the fit
  LKinfo<-  LKrig.setup( x, nlevel=1, alpha=1, NC=30, a.wght=5, lambda=.01)
  obj2<- LKrig( x,y, LKinfo=LKinfo) 

# maximize likelihood over lambda see help( LKrig.MLE) for details 
# this assumes the value of 5 for a.wght.
# In many cases the fit is not very sensitive to the range parameter such as 
# a.wght in this case. 

  MLE.fit<- LKrig.MLE(x,y, LKinfo=LKinfo)
  MLE.fit$out # summary of optimization over lambda.
# fit using MLE for lambda
  obj3<- LKrig( x,y, LKinfo=MLE.fit$LKinfo.MLE)  
  print( obj3)  

# find prediction standard errors at locations based on fixing covariance at MLE's
  out.se<- predict.se( obj3, xnew= x)
# one could evalute the SE on a grid to get the surface of predicted SE's 

##############################################################################################
# Use multiresolution model that approximates an exponential covariance
# with a range of 3 (units of degrees) 
##############################################################################################

  LKinfo<-  Matern2LKrig( range= 3, smoothness=.5, x=x, nlevel=3, NC=10)$LKinfo
# the values of alpha and a.wght that make this work
  LKinfo$alpha
  LKinfo$a.wght
# take a look at the implied covariance function solid=along x dashed=along y
  check.cov<- LKrig.cov.plot( LKinfo)
  matplot( check.cov$d, check.cov$cov, type="l", lty=c(1,2))  
# fit surface
  obj<- LKrig(x,y, LKinfo=LKinfo, lambda=.05)
# summary of fit
  obj
# take a look at the estimate
\dontrun{
  surface(obj)
  points( obj$x)
}
# look at residuals 
  plot( obj$fitted.values, obj$residuals)
#
##############################################################################################
# Search over lambda to find MLE for ozone data with approximate exponential covariance
##############################################################################################

# see help( MLE.search.LKrig) for details 
# recall lamabda controls the choice of measurement error variance (nugget variance)
# see help( Matern2LKrig) for more details about approximate lattice model for exponential
  LKinfo.temp<-  Matern2LKrig( range= 3, smoothness=.5, x=x, nlevel=3, NC=10)$LKinfo
# starting value for lambda optimzation
  LKinfo.temp$lambda<- 1 
  MLE.search<- LKrig.MLE(x,y, LKinfo=LKinfo.temp)
  MLE.ozone.fit<- LKrig( x,y,  LKinfo= MLE.search$LKinfo.MLE)

##############################################################################################
# Including a covariate (linear fixed part in spatial model)
##############################################################################################

  data(COmonthlyMet)
  y.CO<- CO.tmin.MAM.climate
  good<- !is.na( y.CO)
  y.CO<-y.CO[good]
  x.CO<- as.matrix(CO.loc[good,])
  Z.CO<- CO.elev[good]
# single level with large range parameter -- similar to a thin plate spline
#  lambda specified 

# fit with elevation
  obj.CO.elev<- LKrig( x.CO,y.CO,Z=Z.CO, nlevel=1, NC=20,alpha=1, lambda=.005, a.wght=4.1)
# BTW the coefficient for the linear term for elevation  is obj.CO$d.coef[4]
# fitted surface without the elevation term
\dontrun{
   LKinfo<- LKrig.setup( x.CO, nlevel=1, NC=20,alpha=1, a.wght=4.1, lambda=1.0)
  # lambda is the starting vlaue for MLE optimization
  CO.MLE<- LKrig.MLE( x.CO,y.CO,Z=Z.CO, LKinfo=LKinfo)
  obj.CO.elev<- LKrig( x.CO,y.CO,Z=Z.CO, LKinfo= CO.MLE$LKinfo.MLE)
  CO.surface2<- predict.surface( obj.CO.elev, drop.Z=TRUE, nx=50, ny=50)
# pull off CO elevations at same locations on grid as the surface
  data( RMelevation) # a superset of elevations at 4km resolution
  elev.surface<- interp.surface.grid( RMelevation, CO.surface2)
# linear adjustment for elevation 
  CO.surface.elev<- elev.surface$z*obj.CO.elev$d.coef[4]
  CO.full<- CO.surface2$z + CO.surface.elev
# for comparison fit without elevation as a linear covariate:
  CO.MLE2<- LKrig.MLE( x.CO,y.CO, LKinfo=LKinfo)
  obj.CO<- LKrig( x.CO,y.CO, LKinfo= CO.MLE2$LKinfo.MLE)
# surface estimate
  CO.surface<- predict.surface( obj.CO, nx=50, ny=50)
  set.panel( 2,1)
  coltab<- topo.colors(256)
  zr<- range( c( CO.full), na.rm=TRUE)
  image.plot( CO.surface,  col=coltab, zlim =zr)
    US( add=TRUE,lwd=2)
    title( "MAM min temperatures without elevation")
  image.plot( CO.surface$x, CO.surface$y, CO.full, col=coltab, zlim=zr)
    title( "Including elevation")
    US( add=TRUE,lwd=2)
}

#####################################################################################
# for a more elaborate search over the range and lambda to find joint MLEs
# see help( MLE.search.LKrig)
#####################################################################################

#####################################################################################
# A bigger problem: fitting takes about 30 seconds on fast laptop
#####################################################################################
\dontrun{
  data(CO2)
  obj1<- LKrig( CO2$lon.lat,CO2$y,NC=100,nlevel=1, lambda=5, a.wght=5, alpha=1)
# 4600 basis functions 100X46 lattice  (number of basis functions reduced in 
# y direction because of a rectangular domain
  obj1$trA.est # about 1040 effective degrees of freedom 
#
  glist<- list( x= seq( -180,180,,200),y=seq( -80,80,,100) )
  xg<-  make.surface.grid(glist)
  fhat<- predict( obj1, xg)
  fhat <- matrix( fhat,200,100) # convert to image
#Plot data and gap-filled estimate
  set.panel(2,1)
  quilt.plot(CO2$lon.lat,CO2$y,zlim=c(373,381))
  title("Simulated CO2 satellite observations")
  world(add=TRUE,col="magenta")
  image.plot( glist$x, glist$y, fhat,zlim=c(373,381))
  world( add=TRUE, col="magenta")
  title("Gap-filled global predictions")
}  

set.panel()
#####################################################################################
#  Spatially varying alpha example
#####################################################################################
# create a varying alpha


#####################################################################################
#  Comparing LKrig to ordinary Kriging
#####################################################################################

# Here is an illustration of using the fields function mKrig with the LKrig covariance
# to reproduce the computations of LKrig. The difference is that mKrig can
# not take advantage of any sparsity in the precision matrix 
# because its inverse, the covariance matrix, is not sparse.
# This example reinforces the concept that LKrig finds the
# the standard geostatistical estimate but just uses a
# particular covariance function defined via basis functions and 
# the precision matrix and different linear algebra. 
#
# Load ozone data set (AGAIN)
\dontrun{
  data(ozone2)  
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
# Find location that are not 'NA'.
# (LKrig is not set up to handle missing observations.)
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
  a.wght<- 5
  lambda <-  1.5
  obj1<- LKrig( x,y,NC=16,nlevel=1, alpha=1,  lambda=lambda, a.wght=5, NtrA=20,iseed=122)
 
# in both calls iseed is set the same so we can compare 
# Monte Carlo estimates of effective degrees of freedom
  obj1$trA.est
# The covariance "parameters" are all in the list LKinfo
# to create this special list outside of a call to LKrig use
  LKinfo.test <- LKrig.setup( x, NC=16, nlevel=1, alpha=1.0,  a.wght=5)

# this call to mKrig should be identical to the LKrig results
# because it uses the LKrig.cov covariance with all the right parameters.
  obj2<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                      cov.args=list( LKinfo=LKinfo.test), NtrA=20, iseed=122)
# compare the two results this is also an
# example of how tests are automated in fields
# set flag to have tests print results
  test.for.zero.flag<- TRUE
  test.for.zero( obj1$fitted.values, obj2$fitted.values,
                  tag="comparing predicted values LKrig and mKrig")
# compare standard errors. 
  se1<- predict.se.LKrig( obj1)
  se2<- predict.se.mKrig(obj2)
  test.for.zero( se1, se2,
                  tag="comparing standard errors for LKrig and mKrig")
}
#####################################################################################
#  Unnormalized marginal variance for a 6X6 basis on [-1,1]X[-1,1]
#  This is an example of why normalization seems important. 
#####################################################################################

#\dontrun{
# specify covariance without normalization note all we need is the corners of domains
# to setup the info list. 
  LKinfo<- LKrig.setup(cbind( c( -1,1), c(-1,1)), NC=6, nlevel=1,
                                 a.wght=4.5,alpha=1, normalize= FALSE,edge=FALSE)  
# 80X80 grid of points 
  xg<- make.surface.grid( list(x=seq(-1,1,,80), y= seq(-1,1,,80))) 
  look<- LKrig.cov( xg, LKinfo=LKinfo,marginal =TRUE) 
# surface plot of the marginal variances of the process. 
  image.plot( as.surface(xg, look)) 
# basis function centers from the first (and only) level
  xp<- make.surface.grid( LKinfo$grid[[1]])
  points(xp)
#}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

