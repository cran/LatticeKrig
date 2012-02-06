library( LatticeKrig)
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( LatticeKrig)
options( echo=FALSE)

##########################################
  test.for.zero.flag<- 1

LKinfo<- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                     a.wght=c(5,6,7), alpha=1,
                     edge=FALSE, normalize=FALSE)

hold<- LKrig.precision( LKinfo, return.B=TRUE)
hold<- spam2full(hold)

test.for.zero( diag(hold), rep( c(5,6,7), LKinfo$mx*LKinfo$my),
                     tag="diagonal elements of precision 3-levels")
hold2<- LKrig.precision( LKinfo, return.B=TRUE, level.index=2)
hold2<- spam2full( hold2)
number.level<-  LKinfo$mx[2]*LKinfo$my[2]
ind2<-ind1<-  (1:number.level) + LKinfo$offset[2]
test.for.zero( c( hold[ind1, ind2]), c(hold2), tag="just level 2 B matrix")

# now test  t(B)%*%B

hold<- LKrig.precision( LKinfo)
hold<- spam2full(hold)
hold2<- LKrig.precision( LKinfo, level.index=2)
hold2<- spam2full( hold2)
number.level<-  LKinfo$mx[2]*LKinfo$my[2]
ind2<-ind1<-  (1:number.level) + LKinfo$offset[2]
test.for.zero( c( hold[ind1, ind2]), c(hold2), tag="just level 2 Q matrix")

# now everything
  LKinfo<- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4, a.wght=c(5,6,7)
                       , alpha=1)
  hold<- LKrig.precision( LKinfo)
  hold<- spam2full(hold)
  hold2<- LKrig.precision( LKinfo, level.index=2)
  hold2<- spam2full( hold2)
  number.level<-  LKinfo$mx[2]*LKinfo$my[2]
  ind2<-ind1<-  (1:number.level) + LKinfo$offset[2]
  test.for.zero( c( hold[ind1, ind2]), c(hold2), tag= "level 2 Q normalization and edge")

  set.seed(123)
  x1<- cbind( runif( 10), runif(10))
  LKinfo<- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                                                 a.wght=c(5,6,7), alpha=6)
  X<- LKrig.basis(x1, LKinfo)
  X<- spam2full(X)
# check on normalization to unit variance at each level
  Q<- LKrig.precision( LKinfo)
  Q<- spam2full( Q)
  look<- (X)%*% solve( Q) %*%t(X)
  marginal.var<- sum(LKinfo$alpha)
  test.for.zero( diag(look), rep( marginal.var,10),
                tag="normalization to unit variance at each level")
  look2<- LKrig.cov( x1, LKinfo= LKinfo, marginal=TRUE)
  test.for.zero( look2, rep(marginal.var,10) ,
                tag="normalization based on logic in LKrig.cov")
# check full covariance matrix
  look3<- LKrig.cov( x1, x1,LKinfo= LKinfo)
  test.for.zero( look3, look, tag="full covariance from matrix expressions")
# Now test w/o normalization
  LKinfo$normalize<- FALSE
  X<- LKrig.basis(x1, LKinfo)
  X<- spam2full(X)
# check on normalization to unit variance at each level
  Q<- LKrig.precision( LKinfo)
  Q<- spam2full( Q)
  look<- (X)%*% solve( Q) %*%t(X)
  look3<- LKrig.cov( x1, x1,LKinfo= LKinfo)
  test.for.zero( look3, look,
             tag="full covariance from matrix expressions w/o normalization")
#
# check of component covariance matrices
   LKinfo<- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                        a.wght=c(5,6,7), alpha=c(4,2,1))
  set.seed(123)
  x1<- cbind( runif( 10), runif(10))
  x2<- cbind(0,0)
  comp<- matrix( NA,10, 3)   
  for ( l in 1:3){
    grid.info<- LKinfo$grid.info
    grid.info$delta<- LKinfo$delta[l]
    LKinfo.temp<- LKrig.setup( grid.info=grid.info,
                         nlevel=1, a.wght=LKinfo$a.wght[l],
                         alpha=1, edge=TRUE) 
    comp[,l]<- LKrig.cov(x1,x2,LKinfo.temp )
  }
  look1<- comp%*%c( LKinfo$alpha)
  look3<- LKrig.cov( x1, x2,LKinfo= LKinfo)
  test.for.zero( look1, look3, tag="comp normalized cov and LKrig.cov")
#
  cat("Done testing LKrig.precision",fill=TRUE)
  options( echo=FALSE)














