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

LKrig.lnPlike <-
function(Mc,Q,y,lambda,residuals,weights, sigma=NA, rho=NA){
# sanity check on lambda being ratio of sigma^2 and rho
# if sigma and rho are passed.
  if( !is.na(rho) ){
      if( (sigma^2/rho) != lambda){
        stop(" sigma, rho and lambda do not match up")}
  }
#
  n<- length(y) 
  m<- dim(Q)[1]

# find log determinant of reg matrix for use in the log likeihood
  lnDetReg <- 2 * sum(log(diag(Mc)))
# log determinant of precision matrix.
  lnDetQ<-  2* sum( log( diag( chol(Q))))
# now apply a miraculous determinant identity (Sylvester''s theorem)
#  det( I + UV) = det( I + VU)    providing UV is square
# or more generally   
#  det( A + UV) = det(A) det( I + V A^{-1}U)
#  or as we use it
#  ln(det( I + V A^{-1}U)) = ln( det(  A + UV)) - ln( det(A))
#
  lnDetCov<- lnDetReg - lnDetQ + (n - m)* log(lambda) -sum( log(weights))
# finding quadratic form
# this uses a shortcut formula for the quadratic form in terms of
# the residuals
  c.mKrig<- weights*residuals/lambda
  quad.form<-   sum(y* c.mKrig )
# MLE estimate of rho and sigma
# these are derived by assuming Y is  MN(  Td, rho*M )  
  rho.MLE <- quad.form/n
  shat.MLE <- sigma.MLE <- sqrt(lambda * rho.MLE)
# the  log profile likehood with  rho.MLE  and  dhat substituted
# leaving a profile for just lambda.
# note that this is _not_  -2*loglike just the log and
# includes the constants
  lnProfileLike <- (-n/2  - log(2*pi)*(n/2)
                      - (n/2)*log(rho.MLE) - (1/2) * lnDetCov) 
# find log likelihood without profiling if sigma and rho have been passed.
# NOTE: this assumes that  lambda== sigma^2/rho
  if( !is.na(rho) ){
    lnLike<- (-quad.form/(2*rho) - log(2*pi)*(n/2) 
                       - (n/2)*log(rho) -    (1/2) * lnDetCov) }
  else{
    lnLike<- NA}
    
  return( list(lnProfileLike=lnProfileLike,rho.MLE=rho.MLE, shat.MLE=shat.MLE,
               sigma.MLE=shat.MLE,
               sigma=sigma, rho=rho, lnLike=lnLike,               
               quad.form=quad.form, lnDetCov=lnDetCov) )
 
}

