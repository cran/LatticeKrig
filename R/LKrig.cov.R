LKrig.cov <-
function( x1, x2=NULL, LKinfo, C=NA, marginal=FALSE){
  PHI1<- LKrig.basis( x1,LKinfo)
  Q<-  LKrig.precision(LKinfo) 
  Qc <-chol(Q)
# note: construction of lattice basis depends on alpha and a.wght  and already normalizes to
# get unit marginal variances. 
  if( !marginal){
      PHI2<- LKrig.basis( x2,LKinfo)
      if (is.na(C[1])) {
        A <- forwardsolve(Qc, transpose = TRUE, t(PHI2), upper.tri = TRUE)
        A <- backsolve(Qc, A)
      return( PHI1%*%A)}
      else{
        A <- forwardsolve(Qc, transpose = TRUE, t(PHI2)%*%C, upper.tri = TRUE)
        A <- backsolve(Qc, A)
        return(PHI1%*%A)}
  }
    if (marginal){
      if( !is.null(x2)){
      stop("x2 should not be passed to find marginal variance")}
      if(LKinfo$normalize){
        return(rep( 1, nrow(x1)))}
      else{
      Qc<-  chol( LKrig.precision(LKinfo) )
      PHI<- LKrig.basis( x1,LKinfo)
      A <- forwardsolve(Qc, transpose = TRUE, t(PHI), upper.tri = TRUE)
      return( c(colSums( A**2))) }
    }
}

