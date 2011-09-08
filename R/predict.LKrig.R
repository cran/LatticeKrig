predict.LKrig <-
function( object, xnew=NULL,Z=NULL,drop.Z=FALSE,...){
  if( is.null(xnew)){
    xnew<- object$x}
  if( is.null(Z)& object$nZ>0){
    Z<- object$Z} 
  NG<- nrow( xnew)
  if( drop.Z|object$nZ==0){
     temp1<-cbind( rep(1,NG), xnew)%*% object$d.coef[object$ind.drift,]}
   else{
     temp1<- cbind( rep(1,NG), xnew,Z)%*% object$d.coef}
    PHIg<-   LKrig.basis( xnew,object$LKinfo)
    temp2<- PHIg%*%object$c.coef
return( temp1 + temp2)
}

