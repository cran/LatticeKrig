LKrig.sim.conditional<- function( obj, M=1, x.grid= NULL, grid.list=NULL, nx=80, ny=80,...,Z.grid=NULL){
# generate grid if not specified
   if(is.null(x.grid)){
     if( is.null( grid.list)){
       grid.list<- fields.x.to.grid( obj$x, nx=nx,ny=ny)}
    x.grid<- make.surface.grid( grid.list)}    
    ghat<- predict(obj, x.grid, Znew=Z.grid)
# now generate the error surface 
# begin block
  g.conditional.draw<- matrix( NA, ncol=M, nrow=nrow(x.grid))
  set.seed(122) # set seed so same results are produced for help file
  N<- nrow( obj$x)    
  X.full<- rbind( obj$x,x.grid)
# make a copy of basis     
#  PHIg<- LKrig.basis( X.full, obj$LKinfo)   
  for( k in 1:M){
    cat(k, " ")
    g.full <-  ( LKrig.sim( X.full, LKinfo= obj$LKinfo, just.coefficients=FALSE)) * sqrt(obj$rho.MLE)
    g.unconditional.data<- g.full[1:N] 
    g.unconditional.grid<- g.full[ -(1:N)]
    y.synthetic.data<- g.unconditional.data +  obj$sigma.MLE* rnorm(N)/ obj$weights
    obj.fit.synthetic<- LKrig( obj$x, y.synthetic.data,
                              LKinfo=  obj$LKinfo,
                             lambda=obj$lambda, Z=obj$Z, weights=obj$weights,...)
    error<-  g.unconditional.grid -  predict( obj.fit.synthetic,
                                       x.grid, drop.Z= is.null(Z.grid), Znew= Z.grid )
    g.conditional.draw[,k]<- ghat + error
  }
#   
  return( list( x.grid= x.grid, ghat= ghat, g.draw= g.conditional.draw))
}    
