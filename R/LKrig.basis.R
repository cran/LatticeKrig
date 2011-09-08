LKrig.basis <-
function(x1, LKinfo, verbose = FALSE, spam.format=TRUE){
# order of Wendland is hardwired
  Korder<- 2
  grid.info<- LKinfo$grid.info
  nlevel<- LKinfo$nlevel
  overlap<-LKinfo$overlap
  normalize<- LKinfo$normalize
#  
# accumulate matrix column by column in PHI  
  PHI<-NULL
  for( j in 1:nlevel){
# loop over levels the last level might be special ...      
      delta<- LKinfo$delta[j]
      grid.list<- list( x= seq( grid.info$xmin, grid.info$xmax,delta),
                        y= seq( grid.info$ymin, grid.info$ymax,delta))
      centers<- make.surface.grid( grid.list)
      if(verbose){
        print(dim(centers))}
#  set the range of basis functions, they are assumed to be zero outside
#  the radius basis.delta
      basis.delta<- delta*overlap
#      
      PHItemp<-Wendland.basis(x1, centers, basis.delta, max.points = NULL, mean.neighbor = 50,                           just.distance=FALSE)
# accumulate new level of basis function.      
      PHI<-cbind( PHI, PHItemp)
    }
#  normalize to unit marginal variance
    if( normalize){
      if( any(LKinfo$a.wght<=4)){
        stop("Can not normalize with a.wght <= 4")}
      Qc<-  chol( LKrig.precision(LKinfo)  )
      A <- forwardsolve(Qc, transpose = TRUE, t(PHI), upper.tri = TRUE)
      if( nrow(x1)>1){
        wght<- c(colSums( A**2))
        PHI<-diag.spam(1/sqrt(wght))%*%PHI}
      else{
        wght<- sum( A**2)
        PHI@entries<- PHI@entries/ sqrt(wght)}
    }
      
# attach  LKinfo list to the matrix to help identify how the basis functions
# are organized. 
  attr(PHI, which="LKinfo")<-LKinfo
  return(PHI)
}

