# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( LatticeKrig)
options( echo=FALSE)


# tests for computing the determinant and quad form
  test.for.zero.flag<- 1
  alpha<- c(1,.5,.5)
  nlevel<-3
  a.wght<-  c(5,5,10)
  lnDet<- function(A){
  sum( log( eigen( A, symmetric=TRUE)$values))}

  data( ozone2)

  x<-ozone2$lon.lat[1:20,]
  y<- ozone2$y[16,1:20]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
#x<- transformx(x, "range")
  N<- length( y)
  lambda <- .8
 set.seed(123)
  weights<- runif(N)
  W<- diag(weights)
# a micro sized lattice so determinant is not too big or small
  obj<- LKrig( x,y,NC=5, weights= weights, lambda=lambda,nlevel=nlevel,alpha=alpha,a.wght=a.wght,
                              NtrA=5,iseed=122)
  LKinfo<- obj$LKinfo
# now check these formulas as implemented in LatticeKrig
    obj0<- mKrig( x,y, weights= weights, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=LKinfo),
                                 NtrA=20, iseed=122)
 
###### check of formula with weights
  PHI<- spam2full(LKrig.basis( x,LKinfo))
  Q <- spam2full(LKrig.precision(LKinfo))
  M1<- PHI%*%solve( Q)%*%t(PHI) +  lambda*solve( W) 
  B1<- (t(PHI)%*%(W/lambda)%*%PHI + Q)
  B2<- (1/lambda) * ( t(PHI)%*%(W)%*%PHI + lambda*Q)
  B3<-  t(PHI)%*%(W)%*%PHI + lambda*Q
  N2<- nrow(Q)
  lnDet( M1)
#  lnDet( B1) - lnDet(Q) - lnDet( W/lambda)
#  lnDet( B2) - lnDet(Q) - lnDet( W/lambda)
 test.for.zero( lnDet( B3) - lnDet(Q) - sum( log( weights))  + (N-N2)*log(lambda),
                lnDet( M1), tag="Direct formula")
 
 test.for.zero( obj$lnDetCov,  obj0$lnDetCov, tag= "lnDetCov for mKrig and LatticeKrig")
 test.for.zero( obj$quad.form,  obj0$quad.form, tag= "quadratic forms for rho hat")
 test.for.zero(  obj0$lnProfileLike, obj$lnProfileLike,
                                tag="Profile Likelihood concentrated on lambda" )
cat("all done with lnPLike with weights", fill=TRUE)
options( echo=FALSE)
