LKrig.MLE<- function( x,y,..., LKinfo, par.grid=NULL, lambda.profile=TRUE, verbose=FALSE){
  LKrig.args<- c( list( x=x, y=y), list( ...))
  # at this point LKinfo has the correct value for the number of multiresolution levels  
  par.grid<- LKrig.make.par.grid( par.grid=par.grid, LKinfo=LKinfo)
  if( verbose){
    print( par.grid)}
# output matrix
  NG<- length( par.grid$alpha)
  if( verbose){
    print( NG)}
  if(length( par.grid$llambda)!=NG) {
      stop( "llambda values not same length as alpha")}
  if(length( par.grid$a.wght)!=NG) {
      stop( "a.wght values not same length as alpha")}
  out<- matrix( NA, nrow=NG,ncol=6)
  dimnames( out)<- list( NULL,
             c("EffDf", "lnProfLike", "GCV", "sigma.MLE", "rho.MLE", "llambda.MLE"))
#  temporary function used in optimizer  
  temp.fn<- function(x){
    lnLike<- do.call( "LKrig",
               c(LKrig.args, list(LKinfo = LKinfo.temp, lambda=exp(x),
               use.cholesky= MC.save, NtrA=0)))$lnProfileLike
    lnLike}
# first fit to get cholesky symbolic decomposition  
   LKinfo.temp<- LKinfo
   LKinfo.temp$a.wght<- (par.grid$a.wght[[1]])
   LKinfo.temp$alpha<- (par.grid$alpha[[1]])
# save decomp  
   MC.save<-  do.call( "LKrig",
                     c(LKrig.args, list(LKinfo = LKinfo.temp,lambda=1.0,
                     NtrA=0)))$MC
# evaluate parameters but do an optimzation over lambda
  lnProfileLike.max<- -1e20
  for ( k in 1:NG){
    llambda.start<- par.grid$llambda[k]
# if starting value is missing use the optimum from previous fit
# this only really makes sense if other parameter have some sort of continuity from k-1 to k.
    if( is.na( llambda.start) & (k !=1) ){
      llambda.start<-llambda.opt}
 # Note each component of alpha and a.wght is also a list!   
    LKinfo.temp$a.wght<- (par.grid$a.wght[[k]])
    LKinfo.temp$alpha<-  (par.grid$alpha[[k]])
#    if( verbose){
#      print( LKinfo.temp)}
    if( lambda.profile){
      look<- optim( llambda.start, temp.fn,
                   method="BFGS",
                   control=list(fnscale=-1, parscale=.1, ndeps=.01, reltol=1e-6))
      llambda.opt<- look$par
      llambda.start<- llambda.opt
      if( verbose){
        if( k ==1){
        # print heading to list what these are  
          cat("k", "llambda.opt ",   "counts for optim", fill=TRUE)}
        cat( k,  look$par, look$value, look$counts, fill=TRUE)}
      }
    else{
      llambda.opt<- llambda.start}
    obj<- do.call( "LKrig",
               c(LKrig.args, list(LKinfo = LKinfo.temp, lambda=exp(llambda.opt),
               use.cholesky= MC.save, NtrA=20)))
# compare to current largest likelihood and update the LKinfo list if bigger.
    if(obj$lnProfileLike > lnProfileLike.max){
      lnProfileLike.max<- obj$lnProfileLike
      LKinfo.MLE<- obj$LKinfo}
#    
    out[k,1] <- obj$trA.est
    out[k,2] <- obj$lnProfileLike
    out[k,3] <- obj$GCV
    out[k,4] <- obj$sigma.MLE
    out[k,5] <- obj$rho.MLE
    out[k,6]<- llambda.opt
    if( verbose){
      cat("  ", k, "eff.df:", out[k,1], "lnProfile like:", out[k,2],"llambda:", out[k,6],  fill=TRUE)}
  }
  return( list(out=out, par.grid=par.grid, LKinfo=LKinfo,LKinfo.MLE= LKinfo.MLE, call= match.call() ))
  }

