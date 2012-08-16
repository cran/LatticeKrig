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

LKrig.setup <-
function(x=NULL, NC, nlevel, grid.info=NULL,
                       lambda=NA, sigma=NA, rho=NA, 
                       alpha=NA, a.wght=NA, overlap=2.5,
                       normalize=TRUE, normalize.level=NULL,
                       edge=TRUE, rho.object=NULL, RadialBasisFunction="WendlandFunction"){
#
# determines the multiresolution basis function indices.
#
  if( is.null(grid.info)){
    if( is.null(x)){
      stop("need to specify x locations")}
# if center range is missing use the locations
    grid.info<- list( xmin= min( x[,1]), xmax= max( x[,1]),
                      ymin= min( x[,2]), ymax= max( x[,2]))
# spacing for grid
    d1<- grid.info$xmax- grid.info$xmin
    d2<- grid.info$ymax- grid.info$ymin
# actual number of grid points is determined by the spacing delta
# determine delta so that centers are equally
# spaced in both axes and NC is the maximum number of grid points
# along the larger range.     
  grid.info$delta<- max(c(d1,d2))/(NC-1)}
# determine spacing in x and y dimensions.  
  delta<- grid.info$delta   
  grid.list<- list( x= seq( grid.info$xmin, grid.info$xmax,delta),
                    y= seq( grid.info$ymin, grid.info$ymax,delta))
  delta.save<- mx<-my<- rep(NA, nlevel)
# loop through multiresolution levels decreasing delta by factor of 2
# and compute number of grid points. 
  mx[1]<- length( grid.list$x)
  my[1]<- length( grid.list$y)
  delta.save[1]<- delta
 
  grid<- list(grid.list)
# loop through levels doubling the grid resolution  
  if( nlevel >1){
      for( j in 2:nlevel){     
        delta<- delta/2
        delta.save[j]<- delta
        grid.list<- list( x= seq( grid.info$xmin, grid.info$xmax,delta),
                    y= seq( grid.info$ymin, grid.info$ymax,delta))
        grid<- c(grid, list(grid.list))
        mx[j]<- length( grid.list$x)
        my[j]<- length( grid.list$y)      
      }    
    }
  offset<- as.integer( c( 0, cumsum(mx*my)))
# repeat a.wght to fill out for all levels. 
  if( length(a.wght)==1){
    a.wght<- as.list( rep(a.wght[1],nlevel)) }
# coerce a.wght to  list  if it is passed as something else (most likely a vector)
  if( !is.list(a.wght)){
    a.wght<- as.list(a.wght)}
# some checks on a.wght
  a.wght.test<- unlist( a.wght)
  if( any( !is.na(a.wght.test)) ){
     if( any( a.wght.test<4)){              
         stop("a.wght must be >=4 if specified")}
     if( any(a.wght.test==4)&normalize){
         stop("normalize must be FALSE if a.wght ==4")}
  }
# coerce alpha to a list if it is passed as something else
  if( !is.list( alpha)){
    alpha<- as.list( alpha)}
    scalar.alpha <- length(unlist(alpha)) == length( alpha)
# Check some details about scaling the basis functions and how they are
# normalized
  scale.basis<- !is.null(rho.object)
  if( scale.basis & !normalize){
    stop("Can not scale an unnormalized basis")}
  if( is.null(normalize.level)){
     normalize.level= rep( normalize, nlevel)}
# set lambda if sigma and rho are passed. 
  if(is.na(lambda[1])){
    lambda<- sigma^2/ rho}
# 
  out<- list( mx=mx, my=my,nlevel=nlevel,delta= delta.save,m= sum( mx*my),
                 offset=offset,grid.info= grid.info, grid= grid,
                 overlap=overlap, alpha=alpha, a.wght=a.wght,
                 lambda= lambda, sigma=sigma, rho=rho,
                 normalize=normalize,normalize.level=normalize.level, edge=edge,
                 scalar.alpha=scalar.alpha,
                 scale.basis=scale.basis, rho.object=rho.object, RadialBasisFunction=RadialBasisFunction)
  class(out)<-"LKinfo"
  return( out)
}

