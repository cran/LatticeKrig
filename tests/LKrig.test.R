# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( LatticeKrig)
options( echo=FALSE)

##########################################
  test.for.zero.flag<- 1
  data( ozone2)
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]

  N<- length( y)
  a.wght<- 5
  lambda <-  1.5
  obj<- LKrig( x,y,NC=16, lambda=lambda, a.wght=a.wght, alpha=1, nlevel=1, NtrA=20,iseed=122)
  LKinfo<- obj$LKinfo
  K<- LKrig.cov( x,x,LKinfo)
  tempM<-  K
  diag(tempM) <- (lambda) + diag(tempM)
  Mi<- solve( tempM)
  T.matrix<- cbind( rep(1,N),x) 
  d.coef0 <-  solve( t(T.matrix)%*%Mi%*%T.matrix, t(T.matrix)%*%Mi%*%y)
  test.for.zero( obj$d.coef, d.coef0, tag="d from LKrig and by hand")
#### this is c for standard Kriging equations as done in mKrig
  temp2<- chol( tempM)
  c.coef0 <- forwardsolve(temp2, transpose = TRUE,
                        (y- T.matrix%*%d.coef0), upper.tri = TRUE)
  c.coef0 <- backsolve(temp2, c.coef0)
### find these using mKrig (still standard Kriging) and the lattice covariance function:
  obj0<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=LKinfo),
                                 NtrA=20, iseed=122)
  test.for.zero( obj0$c, c.coef0, tag="c from mKrig and by hand" )
# we also know that for standard Kriging
# residuals = lambda* c.coef0
# use this to check the initial LatticeKrig result
 test.for.zero( obj0$fitted.values, obj$fitted.values)
 test.for.zero( lambda*obj0$c, (y-obj$fitted.values),
               tag="c from mKrig and from residuals of LatticeKrig (this is big!)" )

########## redo tests with edge=TRUE
  obj<- LKrig( x,y,NC=16, nlevel=1, alpha=1, lambda=lambda, a.wght=a.wght, NtrA=20,iseed=122, edge=TRUE)
  LKinfo<- obj$LKinfo
  K<- LKrig.cov( x,x,LKinfo)
  tempM<-  K
  diag(tempM) <- (lambda) + diag(tempM)
  Mi<- solve( tempM)
  T.matrix<- cbind( rep(1,N),x) 
  d.coef0 <-  solve( t(T.matrix)%*%Mi%*%T.matrix, t(T.matrix)%*%Mi%*%y)
  test.for.zero( obj$d.coef, d.coef0, tag="d from LKrig and by hand edge=TRUE")
  temp2<- chol( tempM)
#
  c.coef0 <- forwardsolve(temp2, transpose = TRUE,
                        (y- T.matrix%*%d.coef0), upper.tri = TRUE)
  c.coef0 <- backsolve(temp2, c.coef0)
### find these using mKrig (still standard Kriging) and the lattice covariance function:
  obj0<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=LKinfo),
                                 NtrA=20, iseed=122)
  test.for.zero( obj0$c, c.coef0, tag="c from mKrig and by hand edge=TRUE" )
# we also know that for standard Kriging
# residuals = lambda* c.coef0
# use this to check the initial LatticeKrig result
 test.for.zero( obj0$fitted.values, obj$fitted.values)
 test.for.zero( lambda*obj0$c, (y-obj$fitted.values),
               tag="c from mKrig and from residuals of LatticeKrig (this is big!) edge=TRUE" )
######### end tests with edge=TRUE

#
# test more complex covariance model:
#
  alpha<- c(1,.5,.2)
  nlevel<-3
  a.wght<-  c(5,5,10)
  lambda<- .1
  obj<- LKrig( x,y,NC=5, lambda=lambda,
                        nlevel=nlevel, alpha=alpha,a.wght=a.wght, NtrA=20,iseed=122)
  LKinfo<- obj$LKinfo

  obj0<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=LKinfo),
                                 NtrA=20, iseed=122)
  test.for.zero( obj0$fitted.values, obj$fitted.values)
  test.for.zero( obj$d.coef, obj0$d, tag= "d from Lattice Krig and mKrig")
#### check marginal variances
   xg<- make.surface.grid( list(x= seq( -90, -85,, 6), y= seq( 38, 42,,5)) )
   PHI1<- LKrig.basis(x, LKinfo)
   PHI2<- LKrig.basis(xg, LKinfo)                  
   Q<- LKrig.precision( LKinfo)
   Ktest1<- PHI1%*%solve(Q)%*%t(PHI2)
   test.for.zero( Ktest1, LKrig.cov( x,xg, LKinfo=LKinfo))
   Ktest2<- PHI2%*%solve(Q)%*%t(PHI2)
   test.for.zero( diag(Ktest2), LKrig.cov( xg, LKinfo=LKinfo, marginal =TRUE), tag="marginal variance")
#### add some spatially varying alpha's
    N<- LKinfo$mx * LKinfo$my
 
# first a sanity check
    alpha.list<-  list( rep( 1,N[1]), rep(.5,N[2]), rep( .2,N[3]))

    LKinfo2<- LKrig.setup( x, NC=5, nlevel=3, alpha= alpha.list, a.wght=5)
   PHI1<- LKrig.basis(x, LKinfo2)
   PHI2<- LKrig.basis(xg, LKinfo2)                  
   Q<- LKrig.precision( LKinfo2)
   Ktest1<- PHI1%*%solve(Q)%*%t(PHI2)
   test.for.zero( Ktest1, LKrig.cov( x,xg, LKinfo=LKinfo2), tag="spatial alpha cov")
   Ktest2<- PHI2%*%solve(Q)%*%t(PHI2)
   test.for.zero( diag(Ktest2), LKrig.cov( xg, LKinfo=LKinfo2, marginal =TRUE), tag="spatial alpha mar var")
# now varying alphas
   set.seed( 678)
   alpha.list<- list(  runif(N[1]), runif(N[2]), runif(N[3]) )
 LKinfo2<- LKrig.setup( x, NC=5, nlevel=3, alpha= alpha.list, a.wght=5)
   PHI1<- LKrig.basis(x, LKinfo2)
   PHI2<- LKrig.basis(xg, LKinfo2)                  
   Q<- LKrig.precision( LKinfo2)
   Ktest1<- PHI1%*%solve(Q)%*%t(PHI2)
   test.for.zero( Ktest1, LKrig.cov( x,xg, LKinfo=LKinfo2), tag="spatial alpha cov tricksy alpha")
   Ktest2<- PHI2%*%solve(Q)%*%t(PHI2)
   test.for.zero( diag(Ktest2), LKrig.cov( xg, LKinfo=LKinfo2, marginal =TRUE), tag="spatial alpha mar var  tricksy alpha")

#### compare estimates of trace
  test.for.zero( obj$trA.info, obj0$trA.info, tag="Monte Carlo traces")

# evaluate predicted values
  glist<- fields.x.to.grid( x,10, 10)
  xg<-  make.surface.grid(glist)
  grid.info<- obj$LKinfo$grid.info
  LKinfo<- obj$LKinfo
# first "by hand"
  Tmatrix<- cbind( rep(1,nrow(xg)), xg)
  yhat0<- Tmatrix%*%obj0$d +
            LKrig.cov( xg,x, LKinfo)%*%obj0$c
  
  PHIg<- LKrig.basis( xg, LKinfo)
  yhat1<- Tmatrix%*%obj$d.coef + PHIg%*%obj$c.coef
  test.for.zero( yhat1, yhat0, tag="predicted values   by hand and  from standard and DRK")
  test.for.zero( yhat1, predict(obj,xg), tag="predicted values from LatticeKrig and by hand"  )
  test.for.zero( predict(obj,xg), predict(obj0,xg),
                           tag="predicted values LatticeKrig and mKrig")

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
# a micro sized lattice so determinant is not too big or small
  obj<- LKrig( x,y,NC=5, lambda=lambda,nlevel=nlevel,alpha=alpha,a.wght=a.wght,
                              NtrA=5,iseed=122)
  LKinfo<- obj$LKinfo
  grid.info<- LKinfo$grid.info
  PHI<- LKrig.basis( x,LKinfo)
  Q <- LKrig.precision(LKinfo)
# coerce to full matrix
  Q<- spam2full(Q)
  Mtest<- PHI%*% (solve( Q)) %*% t( PHI) + diag(lambda, N)
  temp<- t(PHI)%*%PHI + lambda*Q
  A<- Q*lambda
  B1<-  PHI%*% (solve( A)) %*% t( PHI) + diag(1, N)
  B2<-  t(PHI)%*%PHI + A
# the bullet proof application of identity 
  test.for.zero(lnDet( B1),lnDet( B2)- lnDet(A))
  test.for.zero(
                 lnDet( PHI%*% (solve( Q*lambda)) %*% t( PHI) + diag(1, N)),
                 lnDet( t(PHI)%*%PHI + Q*lambda) - lnDet(Q*lambda) )

# now adjusting for lambda factor 
  test.for.zero( lambda*B1, Mtest)
  test.for.zero(lnDet( Mtest), lnDet(B2) - lnDet(lambda*Q) + N*log(lambda) )
  test.for.zero(lnDet( Mtest), lnDet(B2) - lnDet(Q) + (-LKinfo$m + N)*log(lambda) )

# find log determinant of temp using cholesky factors
# applying det identity
   temp<- t(PHI)%*%PHI + lambda*Q
   chol( temp)-> Mc

   lnDetReg <- 2 * sum(log(diag(Mc)))
   lnDetQ<-  2* sum( log( diag( chol(Q))))
   lnDetCov<- lnDetReg - lnDetQ + (-LKinfo$m + N)*log(lambda)
   test.for.zero( lnDetCov, lnDet( Mtest))
   test.for.zero( obj$lnDetCov, lnDet( Mtest), tag="LatticeKrig and direct test of lnDetCov")
#
###### check of formula with weights
  set.seed(123)
  weights<- runif(N)
  W<- diag(weights)
  lambda<- .5
  PHI<- spam2full(LKrig.basis( x,LKinfo))
  Q <- spam2full(LKrig.precision(LKinfo))
  M1<- PHI%*%solve( Q)%*%t(PHI) +  lambda*solve( W) 

  
  B1<- (t(PHI)%*%(W/lambda)%*%PHI + Q)
  B2<- (1/lambda) * ( t(PHI)%*%(W)%*%PHI + lambda*Q)
  B3<-  t(PHI)%*%(W)%*%PHI + lambda*Q
  N2<- nrow(Q)
  hold<- lnDet( M1)
test.for.zero(   lnDet( B1) - lnDet(Q) - lnDet( W/lambda), hold, tag="Det with weights1")
test.for.zero(   lnDet( B2) - lnDet(Q) - lnDet( W/lambda), hold,  tag="Det with weights1=2")
test.for.zero(  lnDet( B3) - lnDet(Q) - sum( log( weights))  + (N-N2)*log(lambda), hold, tag="Det with weights3")







           
# now check these formulas as implemented in LatticeKrig
 data( ozone2)
  x<-ozone2$lon.lat[1:20,]
  y<- ozone2$y[16,1:20]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
#x<- transformx(x, "range")
  N<- length( y)
  lambda <- .8
# a micro sized lattice so determinant is not too big or small
  obj<- LKrig( x,y,NC=5, lambda=lambda,nlevel=nlevel,alpha=alpha,a.wght=a.wght,
                              NtrA=5,iseed=122)
    obj0<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=LKinfo),
                                 NtrA=20, iseed=122)
 
 test.for.zero( obj$lnDetCov,obj0$lnDetCov, tag= "lnDetCov for mKrig and LatticeKrig")
 test.for.zero( obj$quad.form,  obj0$quad.form, tag= "quadratic forms for rho hat")
 test.for.zero(  obj0$lnProfileLike, obj$lnProfileLike,
                                tag="Profile Likelihood concentrated on lambda" )

# repeat tests for weighted measurement errors.
# recopy data to make reading easier
  data( ozone2)
  x<-ozone2$lon.lat[1:20,]
  y<- ozone2$y[16,1:20]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
#x<- transformx(x, "range")
  N<- length( y)
  alpha<- c(1,.5,.5)
  nlevel<-3
  a.wght<-  c(5,5,10)
  lambda <- 2
  N<- length(y)
  set.seed(243)
  weights<- runif(N)*10
 # weights<- rep( 1, N)
  test.for.zero.flag<- 1
  obj<- LKrig( x,y,weights,NC=15,
                    lambda=lambda,alpha=alpha,nlevel=nlevel, a.wght=a.wght, NtrA=5,iseed=122)
        LKinfo<- obj$LKinfo
# compare mKrig and Krig with weights and LatticeKrig
  obj0<- mKrig( x,y,weights, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=LKinfo),
                                 NtrA=20, iseed=122)
 
  obj1<- Krig( x,y,weights=weights, lambda=lambda,GCV=TRUE, m=2,
               cov.function="LKrig.cov", cov.args=list(LKinfo=LKinfo))
            
 test.for.zero( obj0$fitted.values, obj1$fitted.values)
 test.for.zero( predict(obj0), predict(obj1), tag="predicted  values mKrig/Krig  w/weights")
 test.for.zero( obj0$rhohat, obj1$rhohat,tag="compare rhohat for mKrig and Krig with weights")

############ now tests for LatticeKrig

 test.for.zero( obj$fitted.values, obj0$fitted.values)
 test.for.zero( obj$rho.MLE, obj0$rho.MLE)
 test.for.zero( obj$lnDetCov, obj0$lnDetCov)
############# tests using reuse Mc options
 data( ozone2)
  x<-ozone2$lon.lat[1:20,]
  y<- ozone2$y[16,1:20]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
  N<- length(y)
  set.seed(243)
  weights<- runif(N)*10
#x<- transformx(x, "range")
  N<- length( y)
  alpha<- c(1,.5,.5)
  nlevel<-3
  a.wght<-  c(5,5,10)
  lambda <- .8

 obj<- LKrig( x,y,weights=weights,NC=15, lambda=lambda,alpha=alpha,
                    nlevel=nlevel,a.wght=a.wght, return.cholesky=TRUE)
 obj2<- LKrig( x,y,weights=weights,NC=15, lambda=2*lambda,alpha=alpha,
                    nlevel=nlevel,a.wght=a.wght, use.cholesky=obj$Mc)
 obj3<-  LKrig( x,y,weights=weights,NC=15, lambda=2*lambda,alpha=alpha,
                    nlevel=nlevel,a.wght=a.wght, return.cholesky=TRUE)

 test.for.zero( obj3$c.coef, obj2$c.coef, tag="test of LatticeKrig.coef c")
 test.for.zero( obj3$d.coef, obj2$d.coef, tag="test of LatticeKrig.coef d")

 Q<- LKrig.precision(obj3$LKinfo)
 look2<-LKrig.lnPlike(obj3$Mc,Q,y, 2*lambda,obj3$residuals, weights)
 test.for.zero( look2$lnProfileLike, obj3$lnProfileLike)

# all done!
 cat("Done testing LatticeKrig",fill=TRUE)
 options( echo=FALSE)


# SUPPLEMENT: commented out sanity checks for  weighted/unweighted versions of mKrig and Krig
#hold0<-Krig ( x,y,weights=weights,method="user",GCV=TRUE,lambda=1e-3,
#             cov.function="Exp.simple.cov", cov.args=list( theta=300) )
#hold1<-mKrig(x,y,weights, lambda=1e-3,cov.function="Exp.simple.cov", cov.args=list( theta=300))
#test.for.zero( predict(hold0), predict(hold1))









