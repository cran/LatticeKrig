# fields, Tools for spatial data
# Copyright 2004-2009, Institute for Mathematics Applied to Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"predictSE.LKrig" <- function(object, xnew = object$x, 
    Znew = object$Z, verbose = FALSE, ...) {
    if (is.null(object$Mc)) {
        stop("need to include the sparse cholesky decompostion in LKrig object \r\nin calling LKrig set return.cholesky = TRUE")
    }
    if(is.null(object$fixedFunction)){
      stop("Standard errors not supported with out a fixed component")
    }
    # set some local variables
    NG <- nrow(xnew)
    lambda <- object$lambda
    rho <- object$rho.MLE
    sigma2 <- lambda * rho
    weights <- object$weights
    LKinfo <- object$LKinfo
    # figure out if extra covariates should be included
    # and build fixed effects matrices
    # (ind.drift has the indices of just the spatial drift part --- e.g. the linear polynomial)
    if ( is.null(Znew) & (object$nZ > 0) ) {
        Znew <- object$Z
    }
    T.matrix <-do.call(object$fixedFunction, c(
                                  list( x=object$x, Z=object$Z, distance.type = LKinfo$distance.type),
                                  object$fixedFunctionArgs))
    t0 <- t(do.call(object$fixedFunction, c(
                                  list( x=xnew,     Z=Znew,     distance.type = LKinfo$distance.type),
                                  object$fixedFunctionArgs) ))
 ##       Omega <- object$Omega[object$ind.drift, object$ind.drift]
        Omega <- object$Omega
    #
    k0 <- LKrig.cov(object$x, xnew, LKinfo = LKinfo)
    wS <- diag.spam(sqrt(weights))
    PHI <- LKrig.basis(object$x, LKinfo)
    hold <- LKrig.coef(object$Mc, wS %*% LKrig.basis(object$x, LKinfo),
                     sqrt(weights) * T.matrix, sqrt(weights) * k0, lambda, weights)
    # version of 'c'coefficents for usual Kriging model
    #    c.mKrig<-  weights*(k0 - T.matrix%*%hold$d.coef - PHI%*%hold$c.coef)/ lambda
        d.coef <- hold$d.coef
    # colSums used to find multiple quadratic forms  e.g.  diag(t(x) %*%A%*%x) == colSums( x* (A%*%(x)))
    temp1<-  rho * (colSums(t0*(Omega %*% t0)) - 2*colSums(t0*d.coef))
 # find marginal variances -- trival in the stationary case!
    temp0 <- rho *  (LKrig.cov( xnew, LKinfo=LKinfo, marginal=TRUE) -
   colSums(k0*hold$c.mKrig))
 # Add marginal variance to part from estimate
    temp <- temp0 + temp1
    return(sqrt(temp))
}
