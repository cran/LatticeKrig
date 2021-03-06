# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
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

LKrig.coef <- function(GCholesky, wX, wU, wy, lambda,
                       collapseFixedEffect = FALSE,
                       verbose=FALSE) {
    if (length(lambda) > 1) {
        stop("lambda must be a scalar")
    }
    if( !is.null(wU) ){
        A <- forwardsolve(GCholesky, transpose = TRUE, t(wX) %*% wU, 
                              upper.tri = TRUE)
        A <- backsolve(GCholesky, A)
        A <- t(wU) %*% (wU - wX %*% A)/lambda
#   A is  (T^t M^{-1} T)
        b <- forwardsolve(GCholesky, transpose = TRUE, t(wX) %*% wy, upper.tri = TRUE)
        b <- backsolve(GCholesky, b)
        b <- t(wU) %*% (wy - wX %*% b)/lambda
# b is   (T^t M^{-1} y)
# Save the intermediate matrix   (T^t M^{-1} T) ^{-1}
# this the GLS covariance matrix of estimated coefficients
# should be small for this to be efficient code --  e.g the default is 3X3
        Omega <- solve(A)
# GLS estimates
        d.coef <- Omega %*% b
# combine the different fixed effects estimates across replicates.        
        if(  collapseFixedEffect ){
           d.coefMean<- rowMeans( d.coef)
           dimTemp<- dim ( d.coef)
           d.coef<- matrix( d.coefMean,
                            nrow = dimTemp[1],
                            ncol = dimTemp[2])
        }
        residualFixed<- wy - wU %*% d.coef
   }
   else{
       Omega<- NULL
       d.coef<- NULL
       residualFixed<- wy
   }     
# coefficients of basis functions.
    c.coef <- forwardsolve(GCholesky, transpose = TRUE,
                       t(wX) %*% (residualFixed), upper.tri = TRUE)
# This is the formula from the LKrig article to 
# to evaluate  M^{-1}. The W from this formula is absorbed into weighting of
# the observations and fixed part of the model. ( e.g. wy is  W^{1/2}y )
# Note that older versions have a mistake 
# where the factor (1/lambda) has been omitted. quad.form is not needed for the 
# coefficients but is used in computing the likelihood.
#     
    quad.form<-  (1/lambda) * c( colSums(as.matrix(residualFixed^2))  - 
                                    colSums( as.matrix( c.coef^2) ) )
   c.coef <- backsolve(GCholesky, c.coef)
    if( verbose){
    	cat("d.coef: ", d.coef, fill=TRUE)
    	cat( fill=TRUE)
    	cat("c.coef: ", c.coef, fill=TRUE)
    }
    return( list(c.coef = c.coef, d.coef = d.coef,
                 Omega = Omega, quad.form=quad.form,
                 collapseFixedEffect= collapseFixedEffect )
            )
}

