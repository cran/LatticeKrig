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

LKrig.cov.plot<- function( LKinfo, NP=200){
 grid.info<- LKinfo$grid.info
 ux<- seq( grid.info$xmin, grid.info$xmax,,NP)
 uy<- seq( grid.info$ymin, grid.info$ymax,,NP)
 x1<- cbind( ux,  rep(uy[NP/2], NP) )
 x2<- cbind( ux[NP/2], uy[NP/2] )
 d<- c( rdist(x1,x2) )
# evaluate the covariance from the LKinfo object devised from table
# to approximate the Matern in the X direction
 y<-    c( LKrig.cov( x1,x2, LKinfo) ) 
 x1<- cbind( rep(ux[NP/2], NP), uy)
 d2<- c( rdist(x1,x2) )
# same in Y direction    
 y2<- c( LKrig.cov( x1,x2, LKinfo) )
 return( list(d= cbind(d,d2), cov=cbind( y,y2) ) )
}
