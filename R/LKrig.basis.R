# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu, 
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

LKrig.basis <-
function(x1, LKinfo, verbose = FALSE, spam.format=TRUE){
  if( verbose){
    cat("New version of LKrig.basis to save memory", fill=TRUE)}
# order of Wendland is hardwired
  Korder<- 2
  grid.info<- LKinfo$grid.info
  nlevel<- LKinfo$nlevel
  overlap<-LKinfo$overlap
  normalize<- LKinfo$normalize
  scale.basis<- LKinfo$scale.basis
#  
# accumulate matrix column by column in PHI  
  PHI<-NULL
  for( l in 1:nlevel){
# loop over levels the last level might be special ...      
      delta<- LKinfo$delta[l]
      grid.list<- list( x= seq( grid.info$xmin, grid.info$xmax,delta),
                        y= seq( grid.info$ymin, grid.info$ymax,delta))
      centers<- make.surface.grid( grid.list)
      if(verbose){
        print(dim(centers))}
#  set the range of basis functions, they are assumed to be zero outside
#  the radius basis.delta
      basis.delta<- delta*overlap
#      
      PHItemp<-Radial.basis(x1, centers, basis.delta, max.points = NULL, mean.neighbor = 50,                           just.distance=FALSE, RadialBasisFunction=get(LKinfo$RadialBasisFunction))
       if( normalize){
         if( any((LKinfo$a.wght)[[l]]  <= 4 ) ){
           stop("Can not normalize with a.wght <= 4")}

#  find precision matrix at lewvel alpha without multiplying by 1/alpha         
         tempB <- LKrig.MRF.precision(LKinfo$mx[l], LKinfo$my[l], 
            a.wght = (LKinfo$a.wght)[[l]], edge = LKinfo$edge)
         tempB<- spind2spam(tempB)
         Q.level<- t(tempB)%*%tempB
         wght<- LKrig.quadraticform( Q.level, PHItemp) 
         if( nrow(x1)>1){
           PHItemp<-diag.spam(1/sqrt(wght))%*%PHItemp}
         else{
           PHItemp@entries<- PHItemp@entries/ sqrt(wght)}
        }
# accumulate new level of basis function.      
      PHI<-cbind( PHI, PHItemp)
    }
# include a spatially varying multiplication of process.
# wght are in scale of inverse marginal variance of process
      if( scale.basis){
         wght<- c(predict(LKinfo$rho.object,x1))
         if( nrow(x1)>1){
          PHI<-diag.spam(sqrt(wght))%*%PHI}
      else{
          PHI@entries<- PHI@entries* sqrt(wght)}
      } 
# attach  LKinfo list to the matrix to help identify how the basis functions
# are organized. 
  attr(PHI, which="LKinfo")<-LKinfo
  return(PHI)
}

LKrig.quadraticform<- function( Q, PHI){
#   finds     the quadratic forms    PHI_j^T Q.inverse PHI_j  where PHI_j is the jth row of
#   PHI. 
                        nrow<-  dim( PHI)[1]                   
                        if( nrow>1){
                          BLOCKSIZE<- 2000
                          wght<- rep( NA, nrow)
                          counter<-1 
                          while( counter< nrow){
                            ind<- counter: min((counter+BLOCKSIZE), nrow)
                            A <- forwardsolve(
                                          l=chol(Q),
                                          transpose = TRUE, x=t(PHI[ind,]), upper.tri = TRUE)
                            wght[ind]<- c(colSums( A**2))
                            counter<- counter + BLOCKSIZE}
                        }
#Unfortunately the case when there is one row needs to be handled separately.
                        else{
                           A <- forwardsolve(
                                          l=chol(Q),
                                          transpose = TRUE, x=t(PHI), upper.tri = TRUE)
                          wght<- sum( A**2)}                        
                        return( wght)                                               
}


 
  
