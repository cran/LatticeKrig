LKrig.sim.conditional <- function(LKrigObj, M = 1, x.grid = NULL, 
    grid.list = NULL, nx = 80, ny = 80, ..., Z.grid = NULL, seed=42) {
    # generate grid if not specified
    if (is.null(x.grid)) {
        if (is.null(grid.list)) {
            grid.list <- fields.x.to.grid(LKrigObj$x, nx = nx, ny = ny)
        }
        x.grid <- make.surface.grid(grid.list)
    }
    # NOTE the name x.grid may be misleading because it just needs to a column matrix of
    # locations. It need not follow any regualr pattern.
    # now generate the error surfaces
    # begin block
    # create vector of seeds if needed
    if( length(seed)==1){
        seeds<- seed + ( 0:(M-1))
      }
    #
    g.conditional.draw <-    matrix(NA, ncol = M, nrow = nrow(x.grid))
    d.coef.draw<- matrix(NA, ncol = M, nrow = length( LKrigObj$d.coef) )
    N <- nrow(LKrigObj$x)  
    # complete set of locations to evaluate the field must include the observations too
    PHIGrid<- LKrig.basis(x.grid,LKrigObj$LKinfo)
    PHIObs<- LKrig.basis(LKrigObj$x,LKrigObj$LKinfo)
    # predicted field at grid from the actual data

    spatialPart<- (PHIGrid%*% LKrigObj$c.coef)
    fixedPart<- predictLKrigFixedFunction(LKrigObj, xnew=x.grid, Znew = Z.grid)
    ghat <- fixedPart + spatialPart
    for (k in 1:M) {
        cat(k, " ")
        out<- simConditionalDraw( k,  LKrigObj, ghat, x.grid, Z.grid,  PHIGrid, PHIObs, seeds, ...)
        d.coef.draw[,k] <- out$d.coef
        g.conditional.draw[, k] <- out$g.conditional
    }
    #
    return(list(x.grid = x.grid, ghat = ghat, g.draw = g.conditional.draw,
                           d.coef.draw= d.coef.draw))
}

simConditionalDraw <- function(index=1,  LKrigObj, ghat, x.grid, Z.grid, PHIGrid, PHIObs, seeds= 123, ...){
require(LatticeKrig)
        set.seed( seeds[index] )
# generate process at grid and also on the observation locations.
        simCoefficients<- LKrig.sim(LKinfo = LKrigObj$LKinfo, just.coefficients=TRUE)
        g.unconditional.data <-sqrt(LKrigObj$rho.MLE) *PHIObs%*%simCoefficients 
        g.unconditional.grid <-sqrt(LKrigObj$rho.MLE) *PHIGrid%*%simCoefficients 
        # generate a synthetic data set with fixed part set to zero.
        N<- length( LKrigObj$y)
        y.synthetic.data <- g.unconditional.data + LKrigObj$sigma.MLE * 
            rnorm(N)/LKrigObj$weights
        # use LKrig to find the predictions for the xgrid locations
        # NOTE that LKrig will still estimate the fixed part.
        # and it is important to include this part of estimate
        obj.fit.synthetic <- LKrig(LKrigObj$x, y.synthetic.data, LKinfo = LKrigObj$LKinfo, 
                                   lambda = LKrigObj$lambda, Z = LKrigObj$Z,
                                   weights = LKrigObj$weights,
                                   wPHI = LKrigObj$wPHI, use.cholesky = LKrigObj$Mc, ...)
        d.coef <- obj.fit.synthetic$d.coef
        #
        # predict field
        spatialPart<- (PHIGrid%*% obj.fit.synthetic$c.coef)
        fixedPart<- predictLKrigFixedFunction(obj.fit.synthetic, xnew=x.grid, Znew = Z.grid)
        ghat.synthetic<-  fixedPart + spatialPart
        # add prediction error to the condition mean from the actual data
        g.conditional <- ghat + (g.unconditional.grid -  ghat.synthetic)
        # NOTE: sampling variablity for fixed part is built in too
        # because d.coef are estimated and included in the prediction. 
return(
       list( g.conditional = g.conditional, d.coef = d.coef) )
}
