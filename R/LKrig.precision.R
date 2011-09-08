LKrig.precision <-
function( LKinfo, return.B=FALSE)
{
   mx<- LKinfo$mx
   my<- LKinfo$my
   grid.info<- LKinfo$grid.info
   L<- LKinfo$nlevel
   if(L!= length(my)){
     stop("number of levels and mx and my are not consistent")}
   offset<- LKinfo$offset
   alpha<- LKinfo$alpha
   a.wght<- LKinfo$a.wght 
# ind holds non-zero indices and ra holds the values   
   ind<- NULL
   ra<-NULL
# loop over levels the last level might be special ...   
   for( j in 1:L){
      delta.temp<- LKinfo$delta[j]
# evaluate the H matrix at level j.
# each row of this matrix has an "a"  on diagonal and
# -1  for the nearest neighbot basis functions.
      tempB<-LKrig.MRF.precision(mx[j], my[j],a.wght = a.wght[j], edge=LKinfo$edge)
# accumulate the new block in the growing matrix.
#     the indices that are not zero      
      ind<- rbind( ind, tempB$ind+offset[j])
#     the values of the matrix at these matrix locations      
      ra<- c( ra, sqrt(alpha[j])*tempB$ra)   
    }
# dimensions of the full matrix   
# coerced to integer in LKrig.setup
   da<- c( offset[L+1], offset[L+1])
# convert to spam format:
   temp<-  spind2spam( list( ind= ind, ra=ra, da=da))
   if( return.B){
     return(temp)}
   else{
# find precision matrix Q = t(B)%*%B and return
     return( t(temp)%*%(temp) )}
 }

