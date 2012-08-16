# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


# test of radial basis function based on Wendland
# and using sparse formats
# Important check is of the FORTRAN function dfind2d
# that does pairwise distances among points within a specified range.

library(LatticeKrig)
options( echo=FALSE)
 test.for.zero.flag<- 1
set.seed(123)
x1<-  matrix( runif(40), ncol=2)
x2<-  matrix( runif(30), ncol=2)
n1<- nrow(x1)
n2<- nrow(x2)

look1<-Radial.basis(x1,x2, delta=.7, just.distance=TRUE)
look1<- spam2full(look1)
look2<- rdist( x1,x2)
look2[ look2>.7] <-0
test.for.zero( look1,look2)

# test when range varies among different points
delta<- c( rep(.6,10), rep( .3,n2-10))
look1<- Radial.basis(x1,x2, delta=delta,just.distance=TRUE)
look1<- spam2full(look1)
ind<-matrix( delta, nrow=n1,ncol=n2, byrow=TRUE)
look2<- rdist( x1,x2)
look2[ look2> ind] <- 0
test.for.zero( look1,look2)

look1<-Radial.basis(x1,x2, delta=.5)
look1<- spam2full(look1)
look2<- Wendland2.2(rdist( x1,x2)/.5)
test.for.zero( look1,look2)

delta<- c( rep(.6,10), rep( .3,n2-10))
look1<-Radial.basis(x1,x2, delta=delta)
look1<- spam2full(look1)
dind<-matrix( delta, nrow=n1,ncol=n2, byrow=TRUE)
look2<- Wendland2.2(rdist( x1,x2)/dind)
test.for.zero( look1,look2)

cat( "Done with testing Wendland family", fill=TRUE)
options( echo=TRUE)
