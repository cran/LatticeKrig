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

Matern2LKrig<- function( range, smoothness, x=NULL, xlim=c(-1,1), ylim=c(-1,1), nlevel, NC=10,
                              check=FALSE, NP=500, NCtable=10){
  dname<- paste("LKrig.Ctable", NCtable,nlevel, sep= ".")
# load table   
  data(list=dname)
# make local copy of table  
  tab<- get( dname)
  scale.adjust<- NC/NCtable
  if( scale.adjust < 1) {
    stop("can't have NC less than  NCtable")}
# pull off range from x locations
  if( !is.null(x)){
     xlim<- range( x[,1])
     ylim<- range( x[,2])}  
# standardized range is in (-1,1) convert to this scale. 
   convert.scale<- 2/ ( xlim[2]- xlim[1])
   range.lim<- range( tab$par.grid[,1])* scale.adjust
   sm.lim<- range(tab$par.grid[,2])
   rstand<- range* convert.scale * scale.adjust
   if( (rstand< range.lim[1]) | (rstand >range.lim[2])){
     cat("limits for range parameter:", range.lim/convert.scale, fill=TRUE)
     stop( "range not within table limits")}
   if( (smoothness< sm.lim[1]) | (smoothness > sm.lim[2])){
     cat("limits for smoothness parameter:", sm.lim, fill=TRUE)
     stop( "smoothness not within table limits")}
   kappa<- rep( NA, nlevel)
   alpha <-  rep( NA, nlevel)
   for( k in 1:nlevel){
# NOTE old version of alpha used to ccreate table was the precision matrix weights
#       definition changed to 1/alpha or marginal variance of each component
# it is probably better to interpolate table in reciprocal scale in any case.     
     if( length(tab$par.glist[[2]])==1){
       kappa[k]<-  splint( tab$range, tab$kappa[,k],  rstand )
       alpha[k]<-  1/splint( tab$range, tab$alpha[,k],  rstand )} #changed
     else{
       kappa[k]<-  predict( tab$predict.obj, x= cbind( rstand, smoothness), y=tab$kappa[,k])
       alpha[k]<-  1/predict( tab$predict.obj, x=cbind( rstand, smoothness),
                             y=tab$alpha[,k])} # changed
    }
    a.wght<-  4 + 1/kappa^2  
    LKinfo<-  LKrig.setup( cbind( xlim, ylim), nlevel=nlevel, NC=NC,
                             a.wght= a.wght, alpha= alpha,
                             edge=TRUE,normalize=TRUE)
   obj<- list( LKinfo= LKinfo)
  if( check){
    ux<- seq( xlim[1], xlim[2],, NP)
    uy<- seq( ylim[1], ylim[2],,NP)
    x1<- cbind(ux,  rep(uy[NP/2], NP) )
    x2<- cbind( ux[NP/2], uy[NP/2] )
    obj$x1<- x1
    obj$x2<- x2
    obj$d<- c(rdist(x1,x2))
# evaluate the covariance from the LKinfo object devised from table
# to approximate the Matern in the X direction
    obj$y<-  cbind(Matern( c(obj$d)/range, smoothness= smoothness),
                   LKrig.cov( x1,x2, LKinfo) )
    x1<- cbind(rep(ux[NP/2], NP),uy )
    obj$d2<- c(rdist(x1,x2))
# same in Y direction    
    obj$y2<- cbind(Matern( c(obj$d2)/range, smoothness= smoothness),
                   LKrig.cov( x1,x2, LKinfo) )
  }
  return(obj)
}

