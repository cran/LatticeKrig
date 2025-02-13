# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2024
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

LKrigCovWeightedObs <- function(x1, wX, LKinfo) {
# this is a utility routine to find the covariance between the 
# weighted observation vector and the process at the locations x1. 
#  i.e.  wy = wUd + wXc + error   g= PHI1%*%c
# and  c has covariance matrix   inverse (Q)	
# covariance of (g,y)  is PHI1 %*% inverse(Q) %*% t(wX)
#
	PHI1 <- LKrig.basis(x1, LKinfo)
# sparse precision matrix for the basis coeficients	
	Q <- LKrig.precision(LKinfo)
	Qc <- chol(Q, memory = LKinfo$choleskyMemory)
# note: construction of lattice basis depends on alpha and a.wght
# and normalizes the basis 
			A <- forwardsolve(Qc, transpose = TRUE, t(PHI1), upper.tri = TRUE)
			A <- backsolve(Qc, A)
			return(wX %*% A)
	}
