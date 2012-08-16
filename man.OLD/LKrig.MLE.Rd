% corrections DMH 2/03/2012
\name{LKrig.MLE}
\alias{LKrig.MLE}
\alias{LKrig.make.par.grid}
\alias{which.max.matrix}
\alias{which.max.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simple function to search over covariance parameters for Lattice Krig}
\description{
Given a list of different covariance parameters for the Lattice Krig covariance model 
evaluates the profile likelihood and approximate generalized cross-validation function
at each of the parameter settings. 
}
\usage{
LKrig.MLE( x,y,..., LKinfo, par.grid=NULL, lambda.profile=TRUE, verbose=FALSE)
LKrig.make.par.grid(par.grid=NULL, LKinfo = NULL) 
which.max.matrix(z)
which.max.image(obj)}

\arguments{
  \item{x}{
The spatial locations.
}
  \item{y}{
The observations.
}
  \item{par.grid}{
A list with components  llambda, alpha, a.wght giving the different sets of parameters to 
evaluate -- the kth trial has parameters  par.grid\$llambda[k],  par.grid\$alpha[k,] and
 par.grid\$kappa[k,]. 
 Note that par.grid works with _log_ lambda not lambda. See details below for some other 
features of the par.grid argument.
}
\item{lambda.profile}{If TRUE for each set of parameters in par.grid the value of lambda is found that maximizes the likelihood. 
In this case the llambda value is the starting value for the optimizer. If FALSE the profile likelihood is foudn for the parameters
including llambda. If llammbda[k] is NA then the lambda value found from the k-1 maximization is used as a starting value for the k step. (In the source code this is \code{llambda.opt}.) Of course this only makes sense if the other parameters are ordered so that the results for k-1 make sense as a lambda starting value for k.  }

\item{LKinfo}{
   An LKinfo list that specifies  the latticekrig covariance. Usually this is obtained by a call to LKrig.setup or as the 
 component LKinfo from the object returned by a call to LKrig. 
The search sequentially replaces the alpha and a.wght arguments in this
list by the values in par.grid but leaves everything else the same. If par.grid is not passed the parameter values in
LKinfo are used to evaluate the likelihood. This option is most useful if one has fixed values of alpha and a.wght and the
goal is to maximze the likelihood over lambda. 
}

 \item{obj}{An image list with the usual components x, y, and z.}
  \item{verbose}{
If TRUE prints out intermediate results
}
  \item{z}{A matrix}
  \item{\dots}{
   Any arguments to be passed to LKrig. E.g. \code{x}, \code{y}, \code{Z} a covariate or 
       \code{weights}. }

}
\details{
\code{LKrig.MLE}: 
This is a simple wrapper function to do repeated calls to LKrig to evaluate the
profile likelihood and/or to optimize the likelihood over the lambda parameters. The main point 
is that maximization over the lambda parameter (or equivalently for sigma and rho) is the most important and 
should done before consider variation of pther parameters. 
The covariance/model parameters are alpha, kappa, and log lambda and are separate
matrix or vector components of the \code{par.grid} list. The cleanest version of this function would just require the \code{par.grid} list,
however, to be easier to use there are several options to give partial information and let the function itself create the master parameter list. For example, just a search over lambda should be easy and not require creating \code{par.grid} outside the function. To follow this option one would just give an LKinfo list with an added argument for the starting value for the optimizing over lambda.

The profiling refers to the fact that
the parameters sigma ( measurement error standard deviation) and rho (process marginal variance)
are reduced to  their ratio lambda = sigma**2/ rho  and the MLE for rho and sigma constrained by the ratio are substituted into the likelihood.

In the second example below most of the coding is getting the grid of parameters to search 
in the right form.  It is useful to normalize the alpha parameters to sum to one so that 
the marginal varaince of the process is only parameterized by rho.  To make this easy to implement there is the option to specify the alpha parameters in the form of a mixture 
mode so that the components are positive and add to one. (gamma variable below). If a component \code{gamma} is passed then this is assumed to be in the mixture model form and the alpha weights are computed from this. 

At this point we can not offer any simple empirical rules for choosing the search range for lambda.  Some trial and error is usually required. Typically lambda searchs on an equally spaced log scale, however, are most effective. 

\code{which.max.matrix}{Finds the maximum value in a matrix and returns indices}

\code{which.max.image}{Finds the maximum value in an image matrix and returns indices and the corresponding gird values.}

\code{LKrig.make.par.grid} {This is usually an internal function that converts the list of parameters into an more complex data structure used by LKrig.MLE.}

 }

\value{
  \item{out}{A matrix with columns: effective degrees of freedom, ln Profile likelihood, 
Generalized cross-validation function, MLE sigma, MLE rho}
\item{par.grid }{List of parameters used in search. Some parameters might be filled in from 
            \code{LKinfo} and \code{lambda}.}
\item{LKinfo}{ LKinfo list that was either passed or created}
\item{index.MLE}{Index for case that has largest Likelihood value}
\item{index.GCV}{Index for case that has largest GCV value}
\item{LKinfo.MLE}{LKinfo list at the parameters with largest profile likelihood}
\item{lambda.MLE}{Value of lambda from grid with largest profile likelihood.}
\item{call}{Calling sequence for this function.}

\code{which.max.matrix}{Returns a 2 column  matrix with row and column index of maximum.}

\code{which.max.image}{Returns components x,y,z locating maximum value and component in
giving the row and column of maximum in the image matrix.}

\code{LKrig.make.par.grid}{Returns a list with components, alpha, a.wght. Each component is a list where each 
component of the list is a separate set of parameters. This more general format is useful for the nonstationary
case when the parameters alpha might be a list of nlevel matrices.}

}

\author{Douglas Nychka}

\seealso{
\code{\link{LKrig} }
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# fitting summer precip for  sub region of North America
  data(NorthAmericanRainfall)
# rename for less typing
  x<- NorthAmericanRainfall[,c("longitude", "latitude")]
  y<- log(NorthAmericanRainfall[,c("precip")]) # total precip in 1/10 mm for JJA
# cut down the size of this data set so examples run quickly do not use the
# examples also work with full data sets 
  ind<- x[,1] > -90 & x[,2] < 35 #
  x<- x[ind,]
  y<- y[ind]

# quilt.plot( x,y, nx=200, ny=100); world( add=TRUE, lwd=2, col="grey")

##############################################################################################
# searching  over lambda for single level 
##############################################################################################

# This is a "thin plate spline-like" estimate 
  NG<-5 # make this larger ( e.g. 15) for a better estimate
  par.grid<- list( a.wght= rep( 4.05,NG),alpha= rep(1, NG), llambda=  seq(-8,-2,,NG))
  LKinfo<- LKrig.setup(x,NC=15,nlevel=1)
  lambda.search.results<-LKrig.MLE(
                              x,y,LKinfo=LKinfo, par.grid=par.grid, lambda.profile=FALSE)
  lambda.search.results$out
# profile likelihood
  plot( lambda.search.results$out[,1:2], 
         xlab="effective degrees freedom",
         ylab="ln profile likelihood")
# fit at largest likelihood value:
    lambda.MLE.fit<- LKrig( x,y,
                    LKinfo=lambda.search.results$LKinfo.MLE)
# optimizing  Profile likelihood over lambda using optim
# consider 3 values for a.wght (range parameter)
# in this case the log lambdas passed are the starting values for optim.
  NG<-3
  par.grid<- list( a.wght= c( 4.05,4.1,5),alpha= rep(1, NG), llambda= c(-5,NA,NA))
# NOTE: NAs in llambda mean use the previous MLE for llambda as the current starting value. 
  LKinfo<- LKrig.setup(x,NC=15,nlevel=1)
  lambda.search.results<-LKrig.MLE(
                              x,y,LKinfo=LKinfo, par.grid=par.grid,
                              lambda.profile=TRUE, verbose=TRUE)
  lambda.search.results$out
# note first result a.wght = 4.05 is the optimized result for the grid search given above.

##############################################################################################    
# search over two multi-resolution levels varying the  levels of alpha's
##############################################################################################

# NOTE: search ranges found largely by trial and error to make this example work
#       also the grid is quite coarse ( and NC is small) to be quick as a help file example
  Ndes<- 10  # NOTE: this is set very small just to make example run fast
  set.seed(124)
  par.grid<- list()
# create grid of alphas to sum to 1 use a mixture model parametrization
#  alpha1 = (1/(1 + exp(gamma1)) , alpha2 = exp( gamma1) / ( 1 + exp( gamma1))
# 
  par.grid$gamma<- cbind(runif( Ndes, -3,2), runif( Ndes, -3,2))
  par.grid$a.wght<- matrix( 4.5, nrow=Ndes, ncol=3)
# log lambda grid search values
  par.grid$llambda<- runif( Ndes,-5,-3)  
   LKinfo1<- LKrig.setup( x, NC=5, nlevel=3)
# NOTE: a.wght in call is not used. Also a better search is to profile over llambda

 alpha.search.results<- LKrig.MLE( x,y,LKinfo=LKinfo1, par.grid=par.grid, lambda.profile=FALSE)

##############################################################################################
# Viewing the search results
##############################################################################################

# this scatterplot is good for a quick look because  effective degrees of freedom is a useful summary of fit. 
  plot( alpha.search.results$out[,1:2], 
         xlab="effective degrees freedom",
         ylab="ln profile likelihood")
 


\dontrun{
# a more defensible two level model search 
# with profiling over lambda.
#  This takes a few minutes
  Ndes<- 20 
  nlevel<-2 
  set.seed(123)
  par.grid<- list()
# create grid of alphas to sum to 1 use a mixture model parametrization
#  alpha1 = (1/(1 + exp(gamma1)) , alpha2 = exp( gamma1) / ( 1 + exp( gamma1))
  par.grid$gamma<- runif( Ndes, -3,4)
# values for range (a.wght)
  a.wght<-  4 + 1/exp(seq( 0,4,,Ndes))
  par.grid$a.wght<- cbind( a.wght, a.wght)
# log lambda grid search values
  par.grid$llambda<- rep(-4, Ndes)
  LKinfo1<- LKrig.setup( x, NC=15, nlevel=nlevel)
# NOTE: a.wght in call is not used 
 search.results<- LKrig.MLE( x,y,LKinfo=LKinfo1, par.grid=par.grid,
 lambda.profile=TRUE, verbose=TRUE)
# plotting results
set.panel(1,2)
 plot( search.results$out[,1:2], 
         xlab="effective degrees freedom",
         ylab="ln profile likelihood")
 xtemp<- matrix(NA, ncol=2, nrow=Ndes)
 for( k in 1:Ndes){
   xtemp[k,] <- c( (search.results$par.grid$alpha[[k]])[1], (search.results$par.grid$a.wght[[k]])[1] )}
 quilt.plot( xtemp,search.results$out[,2])
#  fit using Tps
 tps.out<- Tps(  xtemp,search.results$out[,2], lambda=0)
 contour( predict.surface(tps.out), lwd=3,add=TRUE)

set.panel()
}

\dontrun{
# evaluating likelhood at  parameters that approximate an
# exponential covariance function with different ranges. 
# this takes a few minutes
  Nrange<- 10
  range.grid<- seq( 2.5,6,,Nrange)
  nlevel<- 3
  a.wght.grid<- alpha.grid<- matrix( NA, nrow=Nrange, ncol=nlevel)
  for(  j in 1:Nrange){
    LKinfo.approx<- Matern2LKrig( range.grid[j], .5, x, NC=15, nlevel=nlevel)$LKinfo
# the parameters are returned in the list of list format
# need to unlist to save as a matrix row. 
    a.wght.grid[j,]<-  unlist(LKinfo.approx$a.wght)
    alpha.grid[j,]<-  unlist(LKinfo.approx$alpha)} 
  LKinfo.ref<-  Matern2LKrig( mean(range.grid), .5, x, NC=15, nlevel=nlevel)$LKinfo
#  
  par.grid<-list()
  par.grid$alpha<- alpha.grid
  par.grid$a.wght<- a.wght.grid 
# specify first start for llambda but then use previous MLE as the start for new range
  par.grid$llambda<- rep(-3, Nrange)
# evaluate likelihood over these ranges but profile over lambda
  range.search.results<- LKrig.MLE( x,y, par.grid=par.grid,LKinfo=LKinfo.ref)
# same kind of plot as above but more compact code  
# profile surface using effective degrees of freedom and the range.
  plot( range.grid, range.search.results$out[,2], 
         xlab="approx exponential range",
         ylab="ln profile likelihood")
# Try other lambda values
 }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
