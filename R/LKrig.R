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

LKrig <- function(x, y = NULL, weights = rep(1, nrow(x)), 
    Z = NULL, LKinfo = NULL, iseed = 123, NtrA = 20, use.cholesky = NULL, 
    return.cholesky = TRUE, NC, nlevel, a.wght, alpha, nu = NULL, 
    lambda = NA, sigma = NA, rho = NA, rho.object = NULL, overlap = 2.5, 
    normalize = TRUE, edge = TRUE, RadialBasisFunction = "WendlandFunction", 
    V = diag(c(1, 1)), distance.type = "Euclidean", verbose = FALSE) {
    # make sure locations are a matrix and get the rows
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    if (any(duplicated(cat.matrix(x)))) 
        stop("locations are not unique see help(LKrig) ")
    # save seed if random number generation happening outside LKrig
    
    # if LKinfo is missing create it from passed arguments
    if (is.null(LKinfo)) {
        LKinfo <- LKrig.setup(x, NC, nlevel = nlevel, lambda = lambda, 
            sigma = sigma, rho = rho, alpha = alpha, nu = nu, 
            a.wght = a.wght, overlap = overlap, normalize = normalize, 
            edge = edge, RadialBasisFunction = RadialBasisFunction, 
            V = V, distance.type = distance.type, rho.object = rho.object)
    }
    if (!is.na(rho) & !is.na(sigma)) {
        lambda <- sigma^2/rho
    }
    # older code can pass LKinfo but leave out lambda here is a fix for that.
    if (!is.na(lambda)) {
        LKinfo$lambda <- lambda
    }
    # makes sure there are no missing values
    if (!is.null(y)) {
        if (any(is.na(y))) 
            stop("Missing values in y should be removed")
    }
    # make sure covariate is a matrix
    if (!is.null(Z)) {
        Z <- as.matrix(Z)
    }
    lambda = LKinfo$lambda
    if (verbose) {
        cat("lambda", lambda, fill = TRUE)
        print(LKinfo)
    }
    # number of basis functions
    m <- LKinfo$m
    # grid dimensions
    mx <- LKinfo$mx
    my <- LKinfo$my
    #  basis.delta<- overlap*delta
    # weighted observation vector
    wy <- sqrt(weights) * y
    # Spatial drift matrix -- assumed to be linear in coordinates. (m=2)
    # and includes Z covariate
    wT.matrix <- sqrt(weights) * LKrig.fixed.component(x, Z, 
        m = 2, distance.type = LKinfo$distance.type)
    nt <- ncol(wT.matrix)
    nZ <- ifelse(is.null(Z), 0, ncol(Z))
    ind.drift <- c(rep(TRUE, (nt - nZ)), rep(FALSE, nZ))
    # Matrix of sum( N1*N2) basis function (columns) evaluated at the N locations (rows)
    # and multiplied by square root of diagonal weight matrix
    # this can be a large matrix if not encoded in sparse format.
    wS <- diag.spam(sqrt(weights))
    wPHI <- wS %*% LKrig.basis(x, LKinfo, verbose = verbose)
    if (verbose) {
        cat("spam class and dim for wPHI", fill = TRUE)
        print(is.spam(wPHI))
        print(dim(wPHI))
    }
    # square root of precision matrix of the lattice process
    #   solve(t(H)%*%H) is proportional to the covariance matrix of the Markov Random Field
    Q <- LKrig.precision(LKinfo)
    
    # variational matrix used to find coefficients of fitted surface and evaluate
    
    if (verbose) {
        cat("spam class and dim for Q", fill = TRUE)
        print(is.spam(Q))
        print(dim(Q))
    }
    ############################################################################################
    # this is the regularized regression matrix that is the key to the entire algorithm:
    ########################################################################################
    M <- t(wPHI) %*% wPHI + lambda * (Q)
    nzero <- length(M@entries)
    if (verbose) {
        cat("Number of nonzero elements:", nzero, fill = TRUE)
    }
    #
    # S-M-W identity can be used to evaluate the data covariance matrix:
    #   M =  PHI%*% solve(t(H)%*%H) %*% t( PHI) + diag( lambda, N)
    # i.e. because temp is sparse, sparse computations can be used to indirectly solve
    # linear systems based on M
    ############################################################################################
    # find Cholesky square root of this matrix
    ############################################################################################
    #  This is where the heavy lifting happens!  temp is in sparse format so
    #   by the overloading is a sparse cholesky decomposition.
    #  if this function has been coded efficiently this step should dominate
    #  all other computations.
    #  If  a previous sparse cholesky decoposition is passed then the
    #  pattern of sparseness is used for the decoposition.
    #  This can speed the computation as the symbolic decomposition part of the
    #  sparse Cholesky is a nontrivial step. The condition is that
    #  the current 'temp' matrix  has the same sparse pattern as that
    #  which resulted in the factorization  cholesky as 'use.cholesky'
    if (is.null(use.cholesky)) {
        Mc <- chol(M)
    }
    else {
        # reuse a previous decomposition to save computation
        #    # first check that it might be from the same
        #    if( sum( temp@colindices - use.cholesky@colindices)!=0){
        #      stop('use.cholesky not the same sparse pattern as  t(wPHI)%*% wPHI + lambda*(Q) ')}
        Mc <- update.spam.chol.NgPeyton(use.cholesky, M)
    }
    # use Mc to find coefficients of estimate
    out1 <- LKrig.coef(Mc, wPHI, wT.matrix, wy, lambda, weights)
    if (verbose) {
        cat("d.coef", out1$d.coef, fill = TRUE)
    }
    # partially fill object list with some components
    object <- list(x = x, y = y, weights = weights, Z = Z, nZ = nZ, 
        ind.drift = ind.drift, LKinfo = LKinfo, lambda = lambda, 
        sigma = sigma, rho = rho)
    # add in components from coefficient estimates
    object <- c(object, out1)
    #
    fitted.values <- predict.LKrig(object, x, Znew = object$Z)
    # For reference:  fitted.values<- (wT.matrix%*%out1$d.coef + wPHI%*%out1$c.coef)/sqrt(weights)
    residuals <- y - fitted.values
    out2 <- LKrig.lnPlike(Mc, Q, y, lambda, residuals, weights, 
        sigma, rho)
    object <- c(object, out2)
    if (verbose) {
        cat("ln ProfileLike", out2$lnProfileLike, fill = TRUE)
    }
    # estimate trace by Monte Carlo if NtrA greater than zero
    if (NtrA > 0) {
        out3 <- LKrig.traceA(Mc, wPHI, wT.matrix, lambda, weights, 
            NtrA, iseed)
        # find GCV
        out3$GCV = (sum(weights * (residuals)^2)/n)/(1 - out3$trA.est/n)^2
    }
    else {
        out3 <- list(trA.est = NA, trA.SE = NA, GCV = NA)
    }
    object <- c(object, out3)
    #
    if (verbose) {
        print(out3)
    }
    #
    # the output object
    # note the if else switch whether to return the big cholesky decomposition
    object <- c(object, list(fitted.values = fitted.values, residuals = residuals, 
        m = m, lambda.fixed = lambda, nonzero.entries = nzero, 
        spatialdriftorder = 2, nt = nt, eff.df = out3$trA.est, 
        call = match.call()))
    if (return.cholesky) {
        object$Mc <- Mc
    }
    # set the class and return.
    class(object) <- "LKrig"
    return(object)
}

