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

LKrig.basis <- function(x1, LKinfo, verbose = FALSE, 
    spam.format = TRUE) {
    grid.info <- LKinfo$grid.info
    nlevel <- LKinfo$nlevel
    overlap <- LKinfo$overlap
    normalize <- LKinfo$normalize
    scale.basis <- LKinfo$scale.basis
    distance.type <- LKinfo$distance.type
    V <- LKinfo$V
    # We will transform (scale and rotate) x matrix of locations by   x%*% t(solve(V))
    #
    # For the radial basis functions this
    # will imply that the Euclidean distance norm used between two vectors X1 and X2
    # is  t(X1-X2) solve(A) (X1-X2)  with A = (V %*%t(V))
    # Here's the deal on the linear transformation V:
    # It should be equivalent to just transforming the locations outside of LatticeKrig.
    # Where ever the observation locations are used
    # they are first transformed with t(solve(V))).
    # Surprisingly this needs to happen in one place below and in LKrig.setup to
    # determine the grid centers.
    #
    # The RBF centers and the delta scale factor, however, assumed to be
    # in the transformed scale and so a  transformation is not needed.
    # see LKrig.setup for how the centers are determined.
    # The concept is that if LKrig was called with the transformed locations
    # ( x%*%t( solve(V)) instead of x
    # and V set to diag( c(1,1) ) one would get the same results.
    # as passing a non identity V.
    #
    # accumulate matrix column by column in PHI
    PHI <- NULL
    x1Transformed<- x1 %*% t(solve(V))
         if( verbose){
          cat(" Dim x1Transformed",  dim( x1Transformed ), fill=TRUE)
        }
    for (l in 1:nlevel) {
        # Loop over levels and evaluate basis functions in that level.
        # Note that all the center information based on the regualr grids is
        # taken from the LKinfo object
        delta <- LKinfo$delta[l]
        grid.list <- LKinfo$grid[[l]]
        centers <- make.surface.grid(grid.list)
        #  set the range of basis functions, they are assumed to be zero outside
        #  the radius basis.delta
        basis.delta <- delta * overlap
        #     
        PHItemp <- Radial.basis(x1Transformed, centers, 
            basis.delta, max.points = NULL, mean.neighbor = 50, 
            just.distance = FALSE, RadialBasisFunction = get(LKinfo$RadialBasisFunction), 
            distance.type = LKinfo$distance.type)
        if( verbose){
          cat(" Dim PHI level", l, dim( PHItemp), fill=TRUE)
        }

        if (normalize) {
          normtime<-  system.time(
            if( LKinfo$fastNormalization){
               wght <- LKrig.normalize.basis.fast(l, LKinfo, x1Transformed)
             }
            else{
               wght <- LKrig.normalize.basis.slow(l, LKinfo, x1Transformed)
             }
           )
           if( verbose ){
                cat( "Level",l, normtime, fill=TRUE)
              }
# now normalize the basis functions by the weights treat the case with one point separately
          if( nrow( x1)>1){
           PHItemp <- diag.spam( 1/sqrt(wght) ) %*% PHItemp
           }
          else{
              PHItemp@entries <- PHItemp@entries/sqrt(wght)
            }
        }
         

        # accumulate new level of basis function.
        PHI <- cbind(PHI, PHItemp)
    }
    # include a spatially varying multiplication of process.
    # wght are in scale of inverse marginal variance of process
    if (scale.basis) {
        wght <- c(predict(LKinfo$rho.object, x1))
        PHI <- diag.spam(sqrt(wght)) %*% PHI
        }
    # attach  LKinfo list to the matrix to help identify how the basis functions
    # are organized.
    attr(PHI, which = "LKinfo") <- LKinfo
    return(PHI)
}


