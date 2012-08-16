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
LKrig.MRF.precision <-
function(mx,my, a.wght, edge=TRUE) {
  m<- mx*my
  kappa2<- a.wght-4
  da<- as.integer( c(m,m))
  I<- as.integer(rep(1:mx,my))
  J<- as.integer(rep((1:my), rep( mx,my)))
# contents of sparse matrix organize as 5 column matrix
  ra<- cbind(rep( 4+kappa2, m), matrix( -1,m, 4))
#  Note: if  a.wght depends on lattice position
# pass a.wght  as a matrix and use  ra<- as.numeric(c( ra,rep( c(a.wght), 4)))
#  order is  center   top, bottom, left right
#  superset of indices for center and nearest neighbors
#  lattices points on boundaries will need to have some points trimmed
  Bj<- c( I    + (J-1)*mx,
         (I-1) + (J-1)*mx,
         (I+1) + (J-1)*mx,
          I    + (J-2)*mx,
          I    +   (J)*mx )
 # fix up boundaries
     if(edge){
       # the corners: fill each separately for clarity
       # corner.indices<- c( 1, mx, 1 + (my-1)*mx, mx*my)
       # upper left, lower left, upper right lower right
       ra[1,] <-             c(1 + kappa2/4,  NA, -.5, NA, -.5)
       ra[mx,]<-             c(1 + kappa2/4, -.5,  NA, NA, -.5)
       ra[1 + (my-1)*mx,]<-  c(1 + kappa2/4,  NA, -.5, -.5,  NA)
       ra[mx*my,]<-          c(1 + kappa2/4, -.5,  NA, -.5,  NA)
       # edges
       for( j in 2:(my-1)){
       # top row  then bottom row 
         ra[ 1+ (j-1)*mx,]<-      c( 2 + kappa2/2, NA, -1, -.5, -.5)
         ra[ mx + (j-1)*mx,]<-    c( 2 + kappa2/2, -1, NA, -.5, -.5)
        }
        for( i in 2:(mx-1)){
       # left side   then right
          ra[ i , ]<-            c( 2 + kappa2/2, -.5, -.5, NA,  -1)
          ra[ i + (my-1)*mx, ]<- c( 2 + kappa2/2, -.5, -.5, -1,  NA)
        }
     }
# find all cases that are in lattice  
  good<- c( rep( TRUE,m),
           (I-1) > 0,
           (I)   < mx,
           (J-2) >=0,
           (J)   <my )
# remove cases that are beyond the lattice and coerce to integer
# also reshape ra as a vector stacking 5 columns
  Bi<- rep(1:m,5)
  Bi<- as.integer(Bi[good])
  Bj<- as.integer(Bj[good])
  ra<- c(ra)[good]     
# spind format, easier to accumulate columns
# see calling function LKrig.precision
  return( list( ind=cbind(Bi,Bj), ra=ra, da=da))
}

