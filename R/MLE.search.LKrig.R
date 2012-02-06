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

which.max.matrix<- function(z){
  ind<- which.max(z)
  ix<- ind%%nrow(z) -1 
  iy<- ind%/% ncol(z)
  return( cbind( ix,iy))}


which.max.image<- function(obj){
  ind.z<- which.max.matrix( obj$z)
  return( list( x= obj$x[ind.z[,1]], y= obj$y[ind.z[,2]], z=obj$z[ ind.z], ind=ind.z) )  }


MLE.search.LKrig<- function( x,y, par.grid=NULL,NC=NULL,nlevel=1,a.wght=NULL, LKinfo=NULL, llambda=NULL,
                             verbose=TRUE, ...){
  if( is.null(par.grid)){
    par.grid<- list( llambda= llambda)}    
  NG<- length( par.grid$llambda)
# fill out covarinance parameters using the base model in LKinfo if they are not specificed
  if(!is.null(NC) & !is.null(a.wght)){
    LKinfo<- LKrig.setup( x,NC=NC, nlevel=nlevel, a.wght=a.wght)}
# at this point LKinfo has the correct vlaue for the number of mulitresolution levels  
  nlevel<- LKinfo$nlevel
  if( is.null( par.grid$alpha)){
          par.grid$alpha<- matrix( LKinfo$alpha, nrow=NG, ncol=nlevel, byrow=TRUE)}
  if( is.null( par.grid$a.wght)){
          par.grid$a.wght<- matrix( LKinfo$a.wght,  nrow=NG, ncol=nlevel, byrow=TRUE)}
# convert to matrices for easier indexing
  par.grid$alpha<- as.matrix( par.grid$alpha)
  par.grid$a.wght<- as.matrix( par.grid$a.wght)
  par.grid$llambda<-as.matrix( par.grid$llambda)  
# output matrix     
  out<- matrix( NA, nrow=NG,ncol=5)
  dimnames( out)<- list( NULL,
                        c("EffDf", "lnProfLike", "GCV", "sigma.MLE", "rho.MLE"))
  LKinfo.temp<- LKinfo
  for ( k in 1:NG){
    lam.temp<- exp(par.grid$llambda[k])
    LKinfo.temp$alpha<- c(par.grid$alpha[k,])
    LKinfo.temp$a.wght<- c(par.grid$a.wght[k,])
    if( k ==1){
    # fit first model  and also to create/save
    #  parts of cholesky decoposition that can be reused (MC component).
      obj<- LKrig( x,y, LKinfo = LKinfo.temp, lambda=exp(par.grid$llambda[1]), ...) 
      MC.save<- obj$MC}
    else{
      obj <- LKrig( x,y, LKinfo = LKinfo.temp, lambda= lam.temp,...,
                             use.cholesky= MC.save)}
    out[k,1] <- obj$trA.est
    out[k,2] <- obj$lnProfileLike
    out[k,3] <- obj$GCV
    out[k,4] <- obj$shat.MLE
    out[k,5] <- obj$rho.MLE
    if( verbose){
      cat("  ", k, out[k,1:3],  fill=TRUE)}
  } 
  index.MLE<- which.max(out[,2])
  index.GCV<- which.max(out[,3])
# LKinfo list at the MLE values
  LKinfo.temp<- LKinfo
  lam.temp<- exp(par.grid$llambda[index.MLE])
  LKinfo.temp$alpha<- c(par.grid$alpha[k,])
  LKinfo.temp$a.wght<- c(par.grid$a.wght[k,])
  return( list(out=out,par.grid=par.grid, LKinfo=LKinfo,
               index.MLE= index.MLE, index.GCV=index.GCV,
               LKinfo.MLE= LKinfo.temp,  lambda.MLE=lam.temp,
               call=match.call()) )
}
