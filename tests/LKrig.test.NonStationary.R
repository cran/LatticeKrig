# LatticeKrig
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library(LatticeKrig)
#options( echo=FALSE)
 test.for.zero.flag<- 1

# small test dataset
set.seed(123)

x<-  matrix( runif(20), ncol=2)
xg<- make.surface.grid( list( x= seq( 0,1,,15), y= seq( 0,1,,15)))
yg<- xg[,1]*10  + xg[,2]*10 +1
#yg<- rep( 1, nrow(xg))
         
rho.obj<- Tps( xg,yg,lambda=0)
LKinfo1<- LKrig.setup( x, alpha=1, a.wght=5, nlevel=1, NC=4, rho=1,
                                rho.obj= rho.obj, normalize=TRUE)
LKinfo0<- LKrig.setup( x, alpha=1, a.wght=5, nlevel=1, NC=4, rho=1,
                                 normalize=TRUE)

cov0<- LKrig.cov( x,x, LKinfo0)
W0<- LKrig.basis( x, LKinfo0)
Q0<- LKrig.precision( LKinfo0)
cov0B<- W0%*%solve( Q0)%*%t( W0)        
cov0B<- spam2full( cov0B)
test.for.zero( cov0, cov0B, "rho==1 sanity check")

rho1<- predict( rho.obj, x)      
cov1<- LKrig.cov( x,x, LKinfo1)         
cov1A<- diag( c(sqrt( rho1)))%*% cov0B %*% diag( c(sqrt( rho1)))
test.for.zero( cov1, cov1A,"cov weighting by rho.obj")                
Q1<-LKrig.precision( LKinfo1)
W1<- LKrig.basis( x, LKinfo1)         
cov2<-   W1 %*% solve(Q1) %*% t(W1)
test.for.zero( cov1, cov2, tag="rho.obj explicit and LKrig.cov")

         
