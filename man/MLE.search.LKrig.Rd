% corrections DMH 2/03/2012
\name{MLE.search.LKrig}
\alias{MLE.search.LKrig}
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
MLE.search.LKrig(x, y, 
      par.grid=NULL,NC=NULL,nlevel=1,a.wght=NULL, LKinfo=NULL, llambda=NULL,
                             verbose=TRUE, ...) 
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
A list with components  lambda, alpha, kappa giving the different sets of parameters to 
evaluate -- the kth trial has parameters  par.grid\$llambda[k],  par.grid\$alpha[k,] and
 par.grid\$kappa[k,].  
}
\item{NC}{
If LKinfo is omitted, the maximum number of basis functions in the coarsest multi-resolution level
 and in the longest dimesion of the rectangular domain containing the observation locations. 
}
\item{nlevel}{
If LKinfo is omitted, the number of multi-resolution levels.
}
\item{a.wght}{
If LKinfo is omitted, value for a.wght.
}

  \item{LKinfo}{
   An LKinfo list that specifies parameters of the covariance. 
The search sequentially replaces the lambda, alpha and kappa arguments in this
list but leaves everything else the same. 
}
 \item{llambda}{A vector of log lambda values to search. This can be used in place of
  the list \code{par.grid} if this is the only parameter to vary.}
 \item{obj}{An image list with the usual components x, y, and z.}
  \item{verbose}{
If TRUE prints out intermediate results
}
  \item{z}{A matrix}
  \item{\dots}{
   Any other arguments to be passed to LKrig. E.g. \code{Z} a covariate or 
\code{weights}, weights are proportional to the reciprocal measurement error variance.

}
}
\details{
This is a simple wrapper function to do repeated calls to LKrig to evaluate the
profile likelihood. The covariance/model parameters are alpha, kappa, and log lambda and are separate
matrix or vector components of the \code{par.grid} list. The cleanest version of this function would just require the \code{par.grid} list, however, to be easier to use there are several options to give partial information and let the function itself create the master parameter list. For example, just a search over lambda should be easy and not require creating \code{par.grid} outside the function. 
Also, for easy use, the function should not require the LKinfo list that describes the multiresolution to be created for simple choices. As alternative to using LKrig.setup to create this list one can just specify \code{NC} and  \code{a.wght} and if the number of multiresolution levels is different from one \code{nlevel}. 

The profiling refers to the fact that
the parameters sigma ( measurement error standard deviation) and rho (process marginal variance)
are reduced to  their ratio lambda = sigma**2/ rho  and the MLE for rho and sigma constrained by the ratio are substituted into the likelihood.


In the second example below most of the coding is getting the grid of parameters to search 
in the right form. It is useful to normalize the alpha parameters to sum to one so that 
the marginal varaince of the process is only parameterized by rho. 

Note that at this point we can not offer any simple empirical rules for choosing the search range for lambda.  Some trial and error is usually required. Typically lambda searchs on an equally spaced log scale, however, are most effective. 

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

}

\author{Douglas Nychka}

\seealso{
\code{\link{LKrig} }
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# fitting summer precip for North America
  data(NorthAmericanRainfall)
# rename for less typing
  x<- NorthAmericanRainfall[,c("longitude", "latitude")]
  y<- log(NorthAmericanRainfall[,c("precip")]) # total precip in 1/10 mm for JJA
# quilt.plot( x,y, nx=200, ny=100)

##############################################################################################
# searching  over lambda for single level 
##############################################################################################

# This is a "thin plate spline-like" estimate 
  llambda<- seq(-9,-4,,8)
  lambda.search.results<- MLE.search.LKrig(
                              x,y,NC=100, nlevel=1, a.wght=4.001,
                              llambda=llambda)
# profile likelihood
  plot( lambda.search.results$out[,1:2], 
         xlab="effective degrees freedom",
         ylab="ln profile likelihood")
# fit at largest likelihood value:
    lambda.MLE.fit<- LKrig( x,y,lambda= lambda.search.results$lambda.MLE,
                    LKinfo=lambda.search.results$LKinfo.MLE)

##############################################################################################    
# search over two multi-resolution levels varying the  levels of alpha's
##############################################################################################

# NOTE: search ranges found largely by trial and error to make this example work
#       also the grid is quite coarse ( and NC is small) to be quick as a help file example
  Nalpha2<- 5
  Nlambda<- 5
  llambda<- seq(-6,-2,,Nlambda)
  alpha1<- 1.0
  alpha2<- exp(seq( -2, 4,,Nalpha2) )
# searching over both alpha's is the same as varying lambda so just set lambda to 1.0
  par.grid<- list()
  temp.grid<- as.matrix(expand.grid(list( alpha2=alpha2, llambda=llambda)))
# log lambda grid search values
  par.grid$llambda<- temp.grid[,2]
# collapse alphas into matrix
  par.grid$alpha<- cbind( rep( 1, nrow(temp.grid)), temp.grid[,1])
# normalize alphas to sum to one so that marginal variances are constant
# for different alpha2's 
  par.grid$alpha<- par.grid$alpha/ rowSums(par.grid$alpha)
#
  LKinfo<- LKrig.setup( x, NC=30, a.wght= 4.1, nlevel=2, alpha=c(1,.001))
# 
 alpha.search.results<- MLE.search.LKrig(
                        x,y,par.grid=par.grid,LKinfo=LKinfo)

##############################################################################################
# Viewing the search results
##############################################################################################

# this scatterplot is good for quick a look because dof is a useful summary of fit. 
  plot( alpha.search.results$out[,1:2], 
         xlab="effective degrees freedom",
         ylab="ln profile likelihood")
# interpolated ln likelihood surface
  x.design<-cbind( par.grid$llambda , log(par.grid$alpha[,2]/par.grid$alpha[,1]) )
  alpha.fit<- Tps(x.design,alpha.search.results$out[,2], lambda=0)  
  alpha.surface.fit <-predict.surface( alpha.fit, nx=150, ny=150)
  image.plot(  alpha.surface.fit, xlab="ln lambda", ylab=" ln alpha_2 /alpha_1")
  points(x.design, pch=16)  
#  contours close to maximum
  alpha.MLE<- which.max.image( alpha.surface.fit)
  contour( alpha.surface.fit, levels= alpha.MLE$z -qchisq(c(.5,.9,.95,.99), df=4),
              col="magenta", lwd=2, add=TRUE)
  points( alpha.MLE, col="grey50", pch=16 , cex=1.5)
# best point on grid
  Ibest<- alpha.search.results$index.MLE
  y.temp<- log(par.grid$alpha[Ibest,2]/par.grid$alpha[Ibest,1])
  x.temp<- par.grid$llambda[Ibest]
  points( x.temp, y.temp, col="grey50", pch=15, cex=1.5)
#
# effective degrees of freedom at MLE added as a contour
  dof.fit<- Tps( x.design, alpha.search.results$out[,1], lambda=0)
  dof.surface.fit <-predict.surface( dof.fit, nx=150, ny=150)
  dof.MLE<- predict( dof.fit, x= cbind(alpha.MLE$x, alpha.MLE$y) )
  contour(  dof.surface.fit, levels= dof.MLE, col="grey",lty=2,lwd=2, add=TRUE)

##############################################################################################
# Fitting the LKrig surface estiamte with approximate MLEs
# from Tps surface of log Profile likelihood
##############################################################################################

  alpha.par<- c(1, exp( alpha.MLE$y))
  alpha.par<- alpha.par/(sum( alpha.par))
  llambda.MLE<- alpha.MLE$x
  LKinfo.temp<- alpha.search.results$LKinfo.MLE
  LKinfo.temp$alpha<-  alpha.par
# refitting Lattice Krig model with these covariance parameters
  precip.MLE.fit<- LKrig( x,y, LKinfo= LKinfo.temp, lambda= exp(llambda.MLE))
\dontrun{
# compare to closest value from grid search: alpha.search.results$LKinfo.MLE$alpha
  test<-  LKrig( x,y, LKinfo= alpha.search.results$LKinfo.MLE,
                       lambda=  alpha.search.results$lambda.MLE)
# As a sanity check compare:
#  ln Profile likes 
      c( precip.MLE.fit$lnProfileLike, test$lnProfileLike,
         alpha.search.results$out[ alpha.search.results$index.MLE,2],alpha.MLE$z)
#  and effective DFs
      c( precip.MLE.fit$trA.est, test$trA.est,
         alpha.search.results$out[ alpha.search.results$index.MLE,1],dof.MLE)
}
# middle two values should be identical all should be "close"
#
# estimating range parameter and lambda together
# a more elaborate search over exponential like covariances 
# to determine the range and lambda 
#  see help( Matern2LKrig) for more details
\dontrun{
  Nrange<- 10
  Nlambda<- 10
  range.grid<- seq( 8,32,,Nrange)
  llambda<-  seq(-5,-1.5,,Nlambda)
  nlevel<- 3
  a.wght.grid<- alpha.grid<- matrix( NA, nrow=Nrange, ncol=nlevel)
  for(  j in 1:Nrange){
    LKinfo.approx<- Matern2LKrig( range.grid[j], .5, x, NC=15, nlevel=nlevel)$LKinfo
    a.wght.grid[j,]<-  LKinfo.approx$a.wght
    alpha.grid[j,]<-  LKinfo.approx$alpha} 
  LKinfo.ref<-  Matern2LKrig( mean(range.grid), .5, x, NC=15, nlevel=nlevel)$LKinfo
# repeat the parameters to fill out the factorial design   
  par.grid<-list()
  par.grid$alpha<- alpha.grid
  par.grid$a.wght<- a.wght.grid 
  for( k in 2:Nlambda){
    par.grid$alpha<- rbind( par.grid$alpha,  alpha.grid)
    par.grid$a.wght<- rbind( par.grid$a.wght,  a.wght.grid)}
  par.grid$llambda<- rep( llambda, rep( Nrange, Nlambda))
# harmless to add the range component -- it is not used in the search
# but this is a good way to keep track of its values 
  par.grid$range<- rep( range.grid,Nlambda)
#
# the grid search you have been waiting to see
  range.search.results<- MLE.search.LKrig(
                        x,y,par.grid=par.grid,LKinfo=LKinfo.ref)
# same kind of plot as above but more compact code  
# profile surface using effective degrees of freedom and the range.
  range.fit<- Tps( cbind( range.search.results$out[,1] , par.grid$range),
                        range.search.results$out[,2], lambda=0)  
  surface(range.fit, xlab="DOF", ylab="range")  
  points( range.fit$x, pch=16)  

 }












}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
