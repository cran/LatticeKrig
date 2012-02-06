# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
MLE.LKrig<- function(x,y,...,par.grid=NULL, verbose=FALSE){
  N<- length(y)

# capture arguments that also work for a direct call to LKrig
# note that  alpha, a.wght and lambda will need to be modified in the
# actual call as they might be vectors instead of scalars. 
  LKrig.args<- c(list(x=x, y=y),  list(...))
  if( verbose){
    print( names(LKrig.args))}
   if(is.null(par.grid)){
    pg<-as.matrix( expand.grid( LKrig.args$alpha, LKrig.args$a.wght,LKrig.args$lambda))}
  else{
    pg<- par.grid}
  NG<- nrow(pg)
  if(verbose){
    print( dim(pg))}
# find main cholesky decoposition to reuse it in the loop below
# Note: to make the code flow clearer this case will be recomputed in loop.
   LKrig.args$alpha <- pg[1,1]
   LKrig.args$a.wght <- pg[1,2]
   LKrig.args$lambda <- pg[1,3]
   Mc.obj<- do.call("LKrig",c(LKrig.args, list(return.cholesky=TRUE)))$Mc
   if(verbose){
     print(Mc.obj@dimension)}
# initialize arrays for loop  
  lnProfileLike<-rho.MLE<-shat.MLE<-trA<-SEtrA<-GCV<-rep( NA,NG)
  for (  k in 1: NG){
    LKrig.args$alpha <- pg[k,1]
    LKrig.args$a.wght <- pg[k,2]
    LKrig.args$lambda <- pg[k,3]
    if( verbose){
      cat( "alpha, a.wght, lambda", pg[k,],  ":", fill=TRUE)}
# here is the call to LKrig (in do.call format so we can manipulate the
# calling arguments within R. will take advantage of sparse pattern

  out1<- do.call("LKrig",c(LKrig.args, use.cholesky=Mc.obj))
  rho.MLE[k] <- out1$rho.MLE
  shat.MLE[k] <-out1$shat.MLE
  lnProfileLike[k] <-out1$lnProfileLike
  trA[k] <- out1$trA.est
  SEtrA[k]<- sqrt(var(c(out1$trA.info))/ length(out1$trA.info))
  GCV[k]<- out1$GCV
   if(verbose){
     cat("trA", trA[k], fill=TRUE)}
  }
     
  out<-  data.frame(pg, 1/sqrt(pg[,2]-4),trA,SEtrA, shat.MLE, rho.MLE, lnProfileLike, GCV)
  names(out)<-  c("alpha", "a.wght", "lambda","lattice range", "trA", "SEtrA", "sigmaMLE",
                                "rhoMLE","lnProfileLike", "GCV")
  return( out)
}  

