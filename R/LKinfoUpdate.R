
LKinfoUpdate<- function( LKinfo, ...){
   LKinfoCall<- as.list(LKinfo$call)
   argList<- list( ...)
   LKinfoCall[1] <- NULL
 # first replace the calling arguments with the values in LKinfo
  for( argName in names( LKinfoCall) ){
    tempArg<-LKinfo[[ argName]] 
    if( !is.null( tempArg) ){
      LKinfoCall[[ argName]] <- tempArg
    }
   }
 # insert new values for arguments.  
   for( argName in names( argList)){
      LKinfoCall[[ argName]] <- argList[[ argName]]}
# now call setup again to update the LKinfo object with new values.   
   do.call( "LKrig.setup" , LKinfoCall )
 }
