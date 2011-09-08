LKrig.sim <-
function(x1, LKinfo, M=1){
  PHI1<- LKrig.basis(x1,LKinfo)
  Q<-  LKrig.precision(LKinfo)
#  
#  Q is precision matrix of the coefficients -- not of the field
#  last step does the multiplication to go from coefficients to evaluating
#  values at the field
#  Q = t(H)%*%H = inv((Sigma)
#  So   Sigma= Hi%*% t(Hi)  
#  find u= t(Hi) %*% N(0,1)   then cov(u) = t(Hi)%*%Hi
#  Hi is upper triangular
#
# snippet of code to test the algebra ...  
#   x<-seq( 0,1,,20); Sigma<- exp(-rdist( x,x)/2.5); Q<- solve( Sigma)
#   chol(Q)->Mc; H<- Mc ; Hi<- solve(H);
#   test.for.zero( Q, t(H)%*%H); test.for.zero(Sigma, Hi%*%t(Hi))   
#   E<- rnorm(20);  u1<- Hi%*% E ;   u2<-backsolve(Mc,E)
#   test.for.zero(u1,u2)  
# 
  Qc<-chol(Q)
  E<- matrix( rnorm( M*LKinfo$m), nrow=LKinfo$m,ncol= M)
  A <- backsolve(Qc,E)
  return( PHI1%*%A)
}

