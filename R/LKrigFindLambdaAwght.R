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



# ####################################################################################
LKrigFindLambdaAwght <- function(x, y, ...,  LKinfo,
                                 use.cholesky=NULL, 
                                 lowerBoundLogLambda =-16,
                                 upperBoundLogLambda = 4,
                                 lowerBoundOmega = -3,
                                 upperBoundOmega =  .75,
                                 factr=1e7,
                                 pgtol=1e-1,
                                 maxit=15,
                                 verbose=FALSE) {
  #require(stats)
  
  # For rectangle omega = log(kappa) = log(sqrt(Awght-4))
  # but will change with other models. 
  # Parts of the LKrig call that will be fixed.  (except updates to LKinfo)                             
  if( any( attr(LKinfo$a.wght,"isotropic") ) == FALSE  ){
      stop( paste(attr(LKinfo$a.wght,"isotropic"), 
           "findAwght only setup to estimate a single a.wght
                 parameter in the model.")
          )
  }
  LKrigArgs <- c(list(x = x, y = y), list( ...),
                 list( LKinfo = LKinfo,
                         NtrA = 0, 
                  getVarNames = FALSE)
                 )
  
  if( verbose){
    cat( "LKrigFindLambdaAwght: Set of LKrigArgs before first call:", names(LKrigArgs ), fill=TRUE)
  }
  # Set up initial values of Awght and omega
  Awght.init <- as.numeric(LKrigArgs$LKinfo$a.wght[1])
  omega.start <- Awght2Omega( Awght.init, LKinfo)
  
  # Set up initial values of lambda and log lambda
  lambda.start <- LKrigArgs$LKinfo$lambda 
 
  if( is.na(lambda.start) ){ 
    llambda.start <- -1 
  } 
  else{
    llambda.start<- log( lambda.start)
    }
# 
  if ( (llambda.start < lowerBoundLogLambda) || 
       (llambda.start > upperBoundLogLambda) || 
       (is.na(llambda.start))
      ){
    stop("Given lambda value is out of bounds.")
  }  
#
  if(verbose){
    cat("LKrigFindLambdaAwght: llambda.start:",  llambda.start, "a.wght.start:", Awght.init, fill=TRUE)
  }
  a.wghtTemp<- omega2Awght(omega.start, LKinfo)
  lambdaTemp<- exp( llambda.start)
  # initial call to likelihood and also to get symbolic decomposition of 
  # the "M" matrix k-- sparsity pattern does not change as lambda, awght are varied.
  LKrigObject <- do.call("LKrig", c(
    LKrigArgs,
    list( 
      use.cholesky = use.cholesky, 
      return.cholesky = TRUE,
      return.wXandwU = TRUE,
      lambda = lambdaTemp,
      a.wght = a.wghtTemp,
      verbose = FALSE)))
# Update the LKrigArgs with cholesky decomp and wU  
 LKrigArgs$use.cholesky<- LKrigObject$Mc
 LKrigArgs$wU<- LKrigObject$wUb
#
  capture.evaluations <-  rbind( c(lambdaTemp,llambda.start,
                                 a.wghtTemp, omega.start,
                                 LKrigObject$sigma2.MLE.FULL,
                                 LKrigObject$tau.MLE.FULL,
                                 LKrigObject$lnProfileLike.FULL) 
                                )
  if(verbose){
    cat("Capture.evaluations first call", fill=TRUE )
    cat("lambda", "log lambda", "a.wght", "omega",
         "sigma2MLE", "tauMLE", "logProfileLike", fill=TRUE)
    cat( capture.evaluations, fill=TRUE)
  }
  
#####  optimze likelihood over log lambda  and over omega =  log( a.wght -4)/2
    capture.env <- environment()
# last two arguments are specific to this objectinve function     
    result <- try(optim(c(llambda.start, omega.start),
                       LambdaAwghtObjectiveFunction, 
                      lower=c(lowerBoundLogLambda,lowerBoundOmega), 
                      upper=c(upperBoundLogLambda,upperBoundOmega), 
                     method="L-BFGS-B",
#                      method="BFGS",
                      control=list(fnscale = -1,factr=factr,
                                    pgtol=pgtol, maxit=maxit,
                                     ndeps = c(.05,.05)),
                                    
                      LKrigArgs=LKrigArgs,
                      capture.env= capture.env,
                      verbose=verbose
                      ))
    if(verbose){
      cat("Results from optimize:", fill=TRUE)
      print( result )
    }
    evalSummary <- !(class( result)== "try-error")
    llambda.MLE <- result$par[1]
    lambda.MLE<- exp( llambda.MLE)
    omega.MLE<- result$par[2]
    a.wght.MLE<- omega2Awght(omega.MLE,  LKrigArgs$LKinfo )
    LKrigArgs$NtrA <- 20
    
    LKrigObject <- do.call("LKrig", c(LKrigArgs,
                                      list(
                                        lambda = lambda.MLE,
                                        a.wght = a.wght.MLE
                                       )
                                     )
    )
  
  ###### end optimze block    
  # save summary results from this set of parameters.
  # Output to be saved     
  out <- rep(NA, 11)
  names( out) <-  c("EffDf", "lnProfLike", "GCV", "tau.MLE", "sigma2.MLE", 
                    "lambda.MLE", "a.wght.MLE", "lnLike", "functionEval", 
                    "gradientEval", "totalEval")
  out[ 1] <- LKrigObject$trA.est
  out[ 2] <- LKrigObject$lnProfileLike.FULL
  out[ 3] <- LKrigObject$GCV
  out[ 4] <- LKrigObject$tau.MLE.FULL
  out[ 5] <- LKrigObject$sigma2.MLE.FULL
  out[ 6] <- lambda.MLE
  out[ 7] <- a.wght.MLE
  out[ 8] <- LKrigObject$lnLike.FULL
  out[ 9] <- result$counts[1]
  out[10] <- result$counts[2]
  out[11] <- nrow(capture.evaluations )
  
  # Name columns  of likelihood eval. 
  dimnames(capture.evaluations)<- list( NULL,
                c("lambda","logLambda","a.wght","omega",
                  "sigma2.MLE", "tau.MLE","lnProfileLike.FULL"))
  return(list(summary = out,
              LKinfo = LKrigObject$LKinfo,
              llambda.start = llambda.start,
              Awght.start=omega2Awght( omega.start,LKrigArgs$LKinfo),
              lambda.MLE = lambda.MLE,
              a.wght.MLE = a.wght.MLE,
              omega.MLE =omega.MLE,
              llambda.MLE=llambda.MLE,
              lnLike.eval = capture.evaluations,
              call = match.call() )
       )        
}

# Define the objective function 
LambdaAwghtObjectiveFunction<- function(PARS, LKrigArgs, capture.env, verbose=FALSE ) {
  lambdaTemp <- exp( PARS[1] )
  a.wghtTemp <-  omega2Awght( PARS[2], LKrigArgs$LKinfo)
  
  hold <- do.call("LKrig",          
                  c(LKrigArgs, list(
                    lambda = lambdaTemp,
                    a.wght = a.wghtTemp)
                  )
  )[c("sigma2.MLE.FULL","tau.MLE.FULL","lnProfileLike.FULL")] 
  rowForCapture <-c( lambdaTemp,PARS[1],
                     a.wghtTemp,PARS[2],
                     hold$sigma2.MLE.FULL,
                     hold$tau.MLE.FULL,
                     hold$lnProfileLike.FULL 
  )
  if( verbose){
    cat( rowForCapture, fill=TRUE )
  }
  lnProfileLike.FULL<- hold$lnProfileLike.FULL 
  temp.eval <- get("capture.evaluations",
                      envir = capture.env )
  assign("capture.evaluations",rbind(temp.eval, rowForCapture),
                      envir = capture.env)
  return(lnProfileLike.FULL)
}

