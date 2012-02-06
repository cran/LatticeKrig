% corrections DMH 2/03/2012
\name{LatticeKrig}
\alias{LatticeKrig}
\alias{LKrig}
\alias{Matern2LKrig}
\alias{predict.LKrig}
\alias{print.LKrig}
\alias{surface.LKrig}
\title{Spatial prediction and inference using a compactly supported
multi-resolution basis and a lattice model for the basis coefficients.
}
\description{
A variation of Kriging with fixed basis functions that uses the Wendland covariance kernel to 
create a regular set of basis functions on a grid. The coefficients of these
basis functions are modeled as a Markov random field (MRF).  
}
\usage{
LKrig(x,y=NULL, weights = rep(1, nrow(x)),Z=NULL,NC,lambda, LKinfo=NULL,
                        grid.info=NULL, alpha=1.0,a.wght= NULL, beta=NULL, nlevel=1,
                        iseed=123,NtrA=20,
                        use.cholesky=NULL, return.cholesky=FALSE,
                        overlap=2.5,normalize=TRUE, edge=TRUE,
                        rho.object=NULL, 
                        verbose=FALSE)
Matern2LKrig( range, smoothness,x = NULL, xlim=c(-1,1), ylim=c(-1,1), nlevel, NC=10,                              check=FALSE, NP=500, NCtable=10)
\method{predict}{LKrig}( object, xnew=NULL,Z=NULL, drop.Z=FALSE, ...)
\method{surface}{LKrig}( object, ...)
\method{print}{LKrig}( x, digits=4, ...)
}
\arguments{
  
  \item{alpha}{Weights of precision matrices for each level of resolution.}

  \item{a.wght}{The "a" spatial autoregessive parameter for a first order Markov Random field
 (see detials  below) . This controls the spatial dependence and must be greater than or equal to 4. For \code{a.wght} = 4 normalize should be FALSE.
  If there are multiple levels this can be a vector. See Details below.}

  \item{beta}{A parameterization of the MRF dependence  beta = -1/a.wght. This parameter has been diabled -- one should use \code{a.wght}. }

  \item{digits}{Number of digits in printed output.}
   
  \item{drop.Z}{If true the fixed part will only be evaluated at the
spatial drift polynomial part of the fixed model. The contribution from the other
covariates ( e.g elevation) in the fixed component will be omitted.}

  \item{edge}{If TRUE an adjustment is made in precision matrices for edge effects.}

  \item{grid.info}{A list with information to construct the regular grid for the basis function centers.
  The components are: \code{xmin, xmax, ymin,ymax, delta}. If not passed this list is created by
  \code{LKrig.setup} from \code{x} and \code{NC}.}

 \item{iseed}{Random seed used in the Monte Carlo technique
  for approximating the effective degrees of freedom (trace of the
  smoothing matrix) for the GCV criterion. If NA, no seed is set. }

 \item{lambda}{The ratio of the nugget variance (called sigma squared in \code{fields} and \code{LatticeKrig}) to the parameter 
  controlling the marginal variance of the process (called rho in
  \code{fields} and \code{LatticeKrig}).} 

 \item{LKinfo}{A list that is an alternate form of specifying the multiresolution basis. This is usually 
        created by the function \code{LKrig.setup}. If NULL, this list is created and returned as 
a component of the LKrig object. }

   \item{NC}{Maximum number of grid points in one dimension. If the
domain is square there will be NC*NC basis functions/ lattice points in the
coarest level of resolution. If domain is rectangular, the smaller dimension will have less than NC
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
    (with the name \code{Mc}). This option is used with a subsequent call to
     \code{LKrig} with \code{use.cholesky}. }

  \item{rho.object}{A prediction object that defines a spatially varying component for the 
           marginal variance of the process. The object should be such that the code 
   \code{predict(rho.object, x1)} will evalute the function at 2-d locations \code{x1}. See details  below.
} 

  \item{verbose}{If \code{TRUE} print out results as they are computed in loop.}
  
  \item{use.cholesky}{Use the symbolic part of the Cholesky decomposition passed in  this argument.}

  \item{weights}{A vector that is proportional to the reciprocal
  variances of the errors. I.e. errors are assumed to be uncorrelated
  with variances sigma^2/weights. }

  \item{xnew}{Matrix of locations for prediction.}
 
  \item{x}{Spatial locations of observations.}
 
  \item{y}{Spatial observations.}
 
  \item{Z}{Linear covariates to be included in fixed part of the model
  that are distinct from the default first order polynomial in \code{x}.}

  \item{\dots}{Additional arguments to pass to generic methods.}

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
If element in  ith precision/covariance matrix  has a marginal variance of say gamma
 then \code{alpha[i]}gamma will be the marginal variance of the ith component. 
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
For \code{nlevel==1} one needs to specify \code{a.wght} and \code{NC}
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
to take in account both the MRF covaraince and the expansion using the
basis functions. The reader is referred to the function \code{LK.cov}
for an explicit code that computes the implied covariance function for
the process f.  Formally if P is the covariance matrix (the inverse of
Q) for the coefficents then the covariance between the field at two
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
multiplier. This component is specificed by the object
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
operation. lambda is best varied on a log scale but is unfortunately  difficult to interpret independent of the particular set of locations, sample size and covariance. A more useful interpretation of lambda is through the effective degrees of freedom and this value is found by default using a Monte Carlo techinique. The effective degrees of freedom will vary with the dimension of the fixed regression part of the spatial model ( typical 3 = constant + linear terms) and the total number of observations. It can be interpreted as the approximate number of "parameters" needed to represent the spatial prediction surface. For a fixed covariance model the spatial estimate at the observation locations can be 
represented as  f hat = A(lambda) y where A is a matrix that does not depend on observations. 
The effective number of degrees of freedom is the trace of A(lambda) in analogy to the least squares regression "hat" matrix and has an inverse, monotonic relationship to lambda. The Monte Carlo method
uses the fact that if e are iid N(0,1)  E( t(e) A(lambda) e) = tr( A(lambda)). 


\strong{Spatial prediction:} 
Given this model, it is possible to compute the conditional expectation of f
given Y and also evaluate the profile likelihood for lambda, alpha, 
a.wght  and obj.rho (if appropriate). Because both the basis functions and
Q are sparse, it is possible to compute the estimates and likelihood for
large numbers of spatial locations. This model has an advantage over covariance tapering or
compactly supported covariance functions (e.g. fastTps), because the
implied covariance functions can have longer range correlations.

\strong{Descriptions of specific functions:}

\code{LKrig} Find spatial process estimate for fixed covariance specificed by 
 nlevel, alpha, a.wght (or beta), NC (or delta) and lambda.

\code{Matern2LKrig} This is an experimental function to approximate a Matern covariance
function with a multiresolution, LKrig type covariance. Currently only tables for 
the exponential and 3 and 4 levels have been computed and \code{NC} must be 10 or larger. Also this is computed for a limited set of range parameters. An error is returned if the range is larger or smaller than the limits coverd by the tables. 
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
\code{LKrig:} 
An LKrig class object with components for evaluating the estimate at arbitrary locations, 
describing the fit and as an option ( with \code{Mc.return=TRUE}) the Cholesky decomposition to allow for fast updating with new values of lambda, alpha and a.wght. 
The "symbolic" first step in the sparse Cholesky decomposition can also be used to compute the sparse Cholesky decomposition for a different positive definite 
matrix that has the same pattern of zeroes. This option is useful in computing the likelihood
under different covariance parameters. 

The returned component \code{LKinfo} is a list with the information that describes the layout of the multiresolution basis functions and the covariance parameters for the MRF. 
(see \code{LK.basis} as an example.)

\code{Lkrig.cov.plot} Returns components d and cov. Each are two column matrices for evaluating the 
covariance from the center of the domain to the edges along the X and Y directions 
with distances in the "d" component and the covariances in the "cov".

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
##############################################################################################

  obj<- LKrig(x,y,NC=30, lambda=.01, a.wght=5)
# Plot fitted surface.
  surface(obj) # see also predict.surface to just evaluate on a grid
  US(add=TRUE)
# Check effective degrees of freedom in surface.
  obj$trA.est
# NOTE: obj$LKinfo  has the full description of the Lattice Krig covariance model. 

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
  surface(obj)
  points( obj$x)
# look at residuals 
  plot( obj$fitted.values, obj$residuals)
#
##############################################################################################
# Search over lambda to find MLE for ozone data
##############################################################################################

# see help( MLE.search.LKrig) for details 
# recall lamabda controls the choice of measurement error variance (nugget variance)
# (make NG much larger to get more accurate determination of MLE)
# good ranges for lambda presently have to be determined by trial and error.
# these work for this data set and this covariance. 
\dontrun{
  llambda.grid<- seq(-8,4,,8)
  LKinfo<-  Matern2LKrig( range= 3, smoothness=.5, x=x, nlevel=3, NC=10)$LKinfo
  MLE.search<- MLE.search.LKrig(x,y, llambda= llambda.grid, LKinfo=LKinfo)
  set.panel(2,1)
  plot( MLE.search$out[,1],  MLE.search$out[,2], xlab="effective degrees of freedom",
             ylab="ln Profile Likelihood")
  MLE.ozone.fit<- LKrig( x,y, lambda= MLE.search$lambda.MLE, LKinfo= LKinfo)
# default estimate of surface see also  help(predict.surface)
  surface(MLE.ozone.fit)
  points( obj$x)}

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
  obj.CO.elev<- LKrig( x.CO,y.CO,Z=Z.CO,NC=30, lambda=.005, a.wght=4.1)
# BTW the coefficient for the linear term for elevation is obj.CO$d.coef[4]
# fitted surface without the elevation term
\dontrun{
  CO.surface2<- predict.surface( obj.CO.elev, drop.Z=TRUE, nx=50, ny=50)
# pull off CO elevations at same locations on grid as the surface
  data( RMelevation) # a superset of elevations at 4km resolution
  elev.surface<- interp.surface.grid( RMelevation, CO.surface2)
# linear adjustment for elevation 
  CO.surface.elev<- elev.surface$z*obj.CO.elev$d.coef[4]
  CO.full<- CO.surface2$z + CO.surface.elev
# for comparison fit without elevation as a linear covariate:
  obj.CO<- LKrig( x.CO,y.CO,NC=25, lambda=.005, a.wght=4.1)
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
  obj1<- LKrig( CO2$lon.lat,CO2$y,NC=100, lambda=5, a.wght=5)
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
  world(add=TRUE,col="magenta")
  title("Simulated CO2 satellite observations")
  image.plot( glist$x, glist$y, fhat,zlim=c(373,381))
  world( add=TRUE, col="magenta")
  title("Gap-filled global predictions")
}  

set.panel()

#####################################################################################
#  Comparing LKrig to ordinary Kriging
#####################################################################################

# Here is an illustration of using the fields function mKrig with the LKrig covariance
# to reproduce the computations of LKrig. The difference is that mKrig can
# not take advantage of any sparsity in the precision matrix 
# because its inverse, the covariance matrix, is not sparse.
# This example reinforces the concept that LKrig finds the
# the standard geostatistical estiamte but just uses a
# particular covariance function defined via basis functions and 
# the precision matrix and different linear algebra. 
#
  a.wght<- 5
  lambda <-  1.5
  obj1<- LKrig( x,y,NC=16, lambda=lambda, a.wght=5, NtrA=20,iseed=122)
 
# in both calls iseed is set the same so we can compare 
# Monte Carlo estimates of effective degrees of freedom
  obj2<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                      cov.args=list( LKinfo=obj1$LKinfo), NtrA=20, iseed=122)
# The covariance "parameters" are all in the list LKinfo
# to create this outside of a call to LKrig use
  LKinfo.test <- LKrig.setup( x, NC=16,  a.wght=5)

# compare the two results this is also an
# example of how tests are automated in fields
# set flag to have tests print results
  test.for.zero.flag<- TRUE
  test.for.zero( obj1$fitted.values, obj2$fitted.values,
                  tag="comparing predicted values LKrig and mKrig")
  test.for.zero( unlist(LKinfo.test), unlist(obj1$LKinfo),
                  tag="comparing two ways of creating covariance list")

#####################################################################################
#  Unnormalized marginal variance for a 6X6 basis.
#####################################################################################

\dontrun{
  ginfo<- list( xmin=-1, xmax=1, ymin=-1, ymax=1, delta= 2/(6-1)) 
# specify covariance without normalization
  LKinfo<- LKrig.setup(grid.info=ginfo, nlevel=1, a.wght=4.5,alpha=1,
                       normalize= FALSE,edge=FALSE) ] 
# 150X150 grid of points 
  xg<- make.surface.grid( list(x=seq(-1,1,,150), y= seq(-1,1,,150))) 
  look<- LKrig.cov( xg, LKinfo=LKinfo,marginal =TRUE) 
# surface plot of the marginal variances of the process. 
  image.plot( as.surface(xg, look)) 
# basis function centers: 
  xgrid<- seq(ginfo$xmin, ginfo$xmax, ginfo$delta) 
  ygrid<- seq(ginfo$ymin, ginfo$ymax, ginfo$delta) 
  points( expand.grid( xgrid, ygrid)) 
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

