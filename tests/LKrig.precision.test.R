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
                     a.wght=c(5,6,7), alpha=c( 1,1,1),
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
                       , alpha=c(1,1,1))
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
                                                 a.wght=c(5,6,7), alpha=c(6,6,6))
  X<- LKrig.basis(x1, LKinfo)
  X<- spam2full(X)
# check on normalization to unit variance at each level
  Q<- LKrig.precision( LKinfo)
  Q<- spam2full( Q)
  look<- (X)%*% solve( Q) %*%t(X)
  marginal.var<- sum(unlist(LKinfo$alpha))
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
  look1<- comp%*%c(unlist( LKinfo$alpha))
  look3<- LKrig.cov( x1, x2,LKinfo= LKinfo)
  test.for.zero( look1, look3, tag="comp normalized cov and LKrig.cov")
#
# check construction with spatial a.wght
  cat("Now check spatial a.wght and alpha", fill=TRUE)
  a.wght<-  4 + (1:25)*.1
  LKinfo <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=1, NC=5,
                        a.wght=list(a.wght),
                        alpha=1, edge=FALSE)
  look<- LKrig.precision( LKinfo=LKinfo, return.B=TRUE)
  look2<- spam2full( look)
  test.for.zero( diag( look2), a.wght, tag="spatial a.wght 1 level")
# three levels
  a.wght<-  list( 4 + (1:16)*.1, 4+ (1:49)*.1, 4+ (1:169)*.1)
  LKinfo <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                        a.wght=a.wght,
                        alpha=c(1,1,1), edge=FALSE)
  look<- LKrig.precision( LKinfo=LKinfo, return.B=TRUE)
  look2<- spam2full( look)
  test.for.zero( diag( look2), unlist(a.wght), tag="spatial a.wght 3 levels")
#
# edge correction
  a.wght<-  4 + (1:25)*.1
  kappa2 <- matrix(a.wght - 4, 5,5)
  LKinfo <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=1, NC=5,
                        a.wght=list(a.wght),
                        alpha=1, edge=TRUE)
  look<- LKrig.precision( LKinfo=LKinfo, return.B=TRUE)
  look2<- spam2full( look)
  temp<- matrix( diag(look2), 5,5)
  temp2<- matrix( a.wght, 5,5)
# corners
  ind<- rbind( c(1,1), c(5,1), c(1,5), c(5,5))
  test.for.zero( temp[ind],  1+ kappa2[ind]/4, tag="1 level corners")
# edges
  ind<- rbind( cbind(2:4,rep(1,3)),  cbind(rep(1,3),2:4), cbind(2:4,rep(5,3)), cbind(rep(5,3),2:4))
  test.for.zero( temp[ind],  2+ kappa2[ind]/2, tag="1 level edges")
# interior
  ind<- cbind( rep( 2:4,4), rep( 2:4, c(4,4,4)))
  test.for.zero( temp[ind],  4+ kappa2[ind], tag="1 level interior")
# 3 levels
  a.wght<-  list( 4 + (1:16)*.1, 4+ (1:49)*.1, 4+ (1:169)*.1)
  LKinfo <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                        a.wght=a.wght,
                        alpha=c(1,1,1), edge=TRUE)
  look<- LKrig.precision( LKinfo, return.B=TRUE)
  look2<- spam2full( look)
  look3<- look2[ 1:169 +LKinfo$offset[3] , 1:169+ LKinfo$offset[3]]
  temp<- matrix( diag( look3), 13,13)
  kappa2 <- matrix( a.wght[[3]]-4, 13,13)
  ind<- rbind( c(1,1), c(13,1), c(1,13), c(13,13))
  test.for.zero( temp[ind],  1+ kappa2[ind]/4, tag="3rd level corners")
# edges
  ind<- rbind( cbind(2:12,rep(1,11)),  cbind(rep(1,11),2:12), cbind(2:12,rep(13,11)),
              cbind(rep(13,11),2:12))
  test.for.zero( temp[ind],  2+ kappa2[ind]/2, tag="1 level edges")
# interior
  ind<- cbind( rep( 2:12,11), rep( 2:12, rep(11,11) ))
  test.for.zero( temp[ind],  4+ kappa2[ind], tag="1 level interior")
# testing  alpha weighting
# one level
  alpha<-  (1:25)*.1
  LKinfo <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=1, NC=5,
                        a.wght=5,
                        alpha=list(alpha), edge=TRUE)
  look<- LKrig.precision( LKinfo, return.B=TRUE)
  look2<- spam2full(look)
   LKinfo2 <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=1, NC=5,
                        a.wght=5,
                        alpha=1, edge=TRUE)
  look3<-  diag( 1/sqrt(alpha))%*%spam2full(LKrig.precision( LKinfo2, return.B=TRUE))
  test.for.zero( look3, look2, tag="1 level spatial alpha")
  look4<-  spam2full(LKrig.precision( LKinfo))
  test.for.zero( t(look3)%*%look3, look4, tag="1 level spatial alpha Q")
# three levels
  alpha<-  list(  (1:16)*.1, (1:49)*.1, (1:169)*.1)
  LKinfo <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                        a.wght=5,
                        alpha=alpha, edge=TRUE)
  look<- LKrig.precision( LKinfo, return.B=TRUE)
  look2<- spam2full( look)
 
  LKinfo2 <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                        a.wght=5,
                        alpha=c(1,1,1), edge=TRUE)
  look3<- spam2full(LKrig.precision( LKinfo2, return.B=TRUE))
  look3<-  diag( 1/sqrt(unlist(alpha)))%*%look3
  test.for.zero( look3, look2, tag=" 3 levels spatial alpha")
  look4<-  spam2full(LKrig.precision( LKinfo))
  test.for.zero( t(look3)%*%look3, look4, tag="3 levels spatial alpha Q")

  cat("Done testing LKrig.precision",fill=TRUE)
  options( echo=FALSE)














