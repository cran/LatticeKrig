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

LatticeKrig<- function(x, y, Z=NULL, nu=1, nlevel=4, a.wght=4.01, NC=NULL,
                        LKinfo=NULL, na.rm=TRUE, tol=.005,  ...){
  # a crisp wrapper where many default values are exercised.
              ind<- is.na(y)
              if( any(ind)){
                      x<- x[!ind,]
                      y<- y[!ind]
                      warning("NAs removed")
                      if( !is.null(Z)){
                        Z<- as.matrix( Z)[!ind,]
                      }
                    }
              if( is.null(NC)){
              N<- length( y)
              Nbasis<- 4^(nlevel)/ 3
              NCtest<- 2*sqrt( N/(Nbasis))
  # NCtest chosen so that   NCtest^2 * ( 1 + 4 + 16 + 64) ~~ number of basis functions
  # will be about 4*N. 
              NC<- round(max(5, NCtest ))
            }
            if( is.null(LKinfo) ){
              LKinfo<- LKrig.setup( x=x, NC=NC, nu=1, nlevel=4,
                                    a.wght=a.wght,...)
            }    
 # find lambda              
              obj<- LKrigFindLambda( x=x,y=y, Z=Z, LKinfo=LKinfo, tol=tol)
              LKrig( x,y,Z=Z, LKinfo=LKinfo, lambda=obj$lambda.MLE) 
            }
