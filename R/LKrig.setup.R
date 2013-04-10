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

LKrig.setup <- function(x = NULL, NC = NULL, NC.buffer = 5, 
    nlevel, grid.info = NULL, lambda = NA, sigma = NA, rho = NA, 
    alpha = NA, nu = NULL, a.wght = NA, overlap = 2.5, normalize = TRUE, 
    normalize.level = NULL, edge = FALSE, rho.object = NULL, 
    RadialBasisFunction = "WendlandFunction", distance.type = "Euclidean", 
    V = diag(c(1, 1)), verbose = FALSE) {
    #
    # determines the multiresolution basis function indices.
    #
    
    # check distance choices
    if (is.na(match(distance.type, c("Euclidean", "cylinder")))) {
        stop("distance type is not supported (or is misspelled!).")
    }
    #further checks for cylinder case
    if (distance.type == "cylinder") {
        if (V[1, 1] != 1) {
            stop("can not scale the angular coordinate (x[,1]),\nassumed to be in degrees")
        }
        if ((V[2, 1] != 0) | (V[1, 2] != 0)) {
            stop("can have off diagonal elements in V")
        }
    }
    #
    if (is.null(grid.info)) {
        # if center grid information is missing create grid based on the locations
        # find range of scaled locations
        if (is.null(x)) {
            stop("need to specify x locations")
        }
        if (is.null(NC)) {
            stop("need to specify NC for grid size")
        }
        range.x <- apply(as.matrix(x) %*% t(solve(V)), 2, "range")
        if (verbose) {
            cat("ranges of transformed variables", range.x, fill = TRUE)
        }
        grid.info <- list(xmin = range.x[1, 1], xmax = range.x[2, 
            1], ymin = range.x[1, 2], ymax = range.x[2, 2])
        d1 <- grid.info$xmax - grid.info$xmin
        d2 <- grid.info$ymax - grid.info$ymin
        grid.info$delta <- ifelse(distance.type == "cylinder", 
            d1/(NC - 1), max(c(d1, d2))/(NC - 1))
    }
    
    # check that cylinder geometry has the right delta with x i.e.divides range evenly.
    if (distance.type == "cylinder") {
        d1 <- grid.info$xmax - grid.info$xmin
        NC.test <- d1/grid.info$delta
        if (round(NC.test - round(NC.test)) > 1e-08) {
            stop("delta spacing in x dimension must be even for\ncylinder geometry")
        }
    }
    # find the grid centers
    out.grid <- LKrig.make.centers(grid.info, nlevel, NC.buffer, 
        distance.type)
    #
    # if a.wght is not a list with nlevel components then
    # assume this is a scalar or vector of values for center value.
    if (!is.list(a.wght)) {
        # some checks on a.wght
        # coerce a.wght to list if it is passed as something else (most likely a vector)
        if (nlevel == 1) {
            a.wght <- list(a.wght)
        }
        else {
            # repeat a.wght to fill out for all levels.
            if (length(a.wght) == 1) {
                a.wght <- rep(a.wght, nlevel)
            }
            a.wght <- as.list(c(a.wght))
        }
    }
    #################################
    # check length of a.wght list
    if (length(a.wght) != nlevel) {
        stop("length of a.wght list differs than of nlevel")
    }
    #
    # now figure out if the model is stationary
    # i.e. a.wght pattern is to be  repeated for each node
    # this is the usual case
    # if not stationary a.wght should lists of arrays that
    # give values for each node separately
    stationary <- is.null(dim(a.wght[[1]]))
    # simple check on sizes of arrays
    if (stationary) {
        N.a.wght <- length(a.wght[[1]])
        # allowed lengths for a.wght are just the center 1 values
        # or 9 values for center,  first, and second order neighbors
        if (is.na(match(N.a.wght, c(1, 9)))) {
            stop("a.wght needs to be of length 1 or 9")
        }
        for (k in 1:length(a.wght)) {
            if (length(length(a.wght[[k]])) != N.a.wght) {
                stop(paste("a.wght level", k, " should have length", 
                  N.a.wght))
            }
        }
    }
    else {
        for (k in 1:length(a.wght)) {
            dim.a.wght <- dim(a.wght[[k]])
            if ((dim.a.wght[1] != out.grid$mx[k]) | (dim.a.wght[2] != 
                out.grid$my[k])) {
                stop(paste("a.wght array at level", k, " has wrong first two dimensions", 
                  ))
            }
        }
    }
    ############################################# end a.wght checks
    if (!is.null(nu)) {
        alpha <- exp(-2 * (1:nlevel) * nu)
        alpha <- alpha/sum(alpha)
    }
    # coerce alpha to a list if it is passed as something else
    if (!is.list(alpha)) {
        alpha <- as.list(alpha)
    }
    scalar.alpha <- length(unlist(alpha)) == length(alpha)
    # Check some details about scaling the basis functions and how they are
    # normalized
    scale.basis <- !is.null(rho.object)
    if (scale.basis & !normalize) {
        stop("Can not scale an unnormalized basis")
    }
    if (is.null(normalize.level)) {
        normalize.level = rep(normalize, nlevel)
    }
    # set lambda if sigma and rho are passed.
    if (is.na(lambda[1])) {
        lambda <- sigma^2/rho
    }
    #
    out <- list(nlevel = nlevel, grid.info = grid.info, NC = NC, 
        NC.buffer = NC.buffer, delta = out.grid$delta.save, mx = out.grid$mx, 
        my = out.grid$my, m = out.grid$m, m.domain = out.grid$m.domain, 
        NC.buffer.x = out.grid$NC.buffer.x, NC.buffer.y = out.grid$NC.buffer.y, 
        offset = out.grid$offset, grid = out.grid$grid, overlap = overlap, 
        alpha = alpha, nu = nu, a.wght = a.wght, stationary = stationary, 
        lambda = lambda, sigma = sigma, rho = rho, normalize = normalize, 
        normalize.level = normalize.level, edge = edge, scalar.alpha = scalar.alpha, 
        scale.basis = scale.basis, rho.object = rho.object, RadialBasisFunction = RadialBasisFunction, 
        distance.type = distance.type, V = V)
    class(out) <- "LKinfo"
    return(out)
}

