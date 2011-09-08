LKrig <-
function( x,y=NULL, weights = rep(1, nrow(x)),Z=NULL,NC,lambda, LKinfo=NULL,
                        grid.info=NULL, alpha=1.0,a.wght= 5, beta=NULL, nlevel=1,
                        iseed=123,NtrA=20,
                        use.cholesky=NULL, return.cholesky=FALSE,
                        overlap=2.5, normalize=TRUE,edge=FALSE,
                        verbose=FALSE){
# make sure locations are a matrix and get the rows  
  x<- as.matrix(x)
  n<- nrow(x)
  if (any(duplicated(cat.matrix(x)))) 
        stop("locations are not unique see help(LKrig) ")
# makes sure there are no missing values
  if(!is.null(y)){
    if (any(is.na(y))) 
      stop("Missing values in y should be removed")}
# make sure covariate is a matrix  
   if(!is.null(Z)){
      Z<- as.matrix(Z)}
# the following function creates the master list LKinfo
# if it has not been passed
# this list describes the multi-resolution covariance model.
  if( is.null(LKinfo)){
    LKinfo<- LKrig.setup(x,NC, grid.info, nlevel=nlevel,  alpha=alpha, a.wght=a.wght, beta=beta,
                            overlap=overlap, normalize=normalize, edge=edge)}
  if(verbose){
    print(LKinfo)}
# number of basis functions   
   m<-  LKinfo$m
# grid dimensions
   mx<- LKinfo$mx
   my<- LKinfo$my
#  basis.delta<- overlap*delta
# weighted observation vector  
  wy<- sqrt(weights)*y
# Spatial drift matrix -- assumed to be linear in coordinates. (m=2)
# and includes Z covariate if not NULL 
  wT.matrix<- sqrt(weights)*cbind( rep(1, n), x,Z)
  nt<- ncol(wT.matrix)
  nZ<- ifelse (is.null(Z), 0, ncol(Z))   
  ind.drift <- c(rep(TRUE, (nt-nZ)), rep( FALSE, nZ))
# Matrix of sum( N1*N2) basis function (columns) evaluated at the N locations (rows)
# and multiplied by square root of diagonal weight matrix  
# this can be a large matrix if not encoded in sparse format.
  wS<-  diag.spam(sqrt( weights))
 
  wPHI<-wS%*%LKrig.basis(x,LKinfo)
  if(verbose){
    cat("spam class and dim for wPHI", fill=TRUE)
    print( is.spam(wPHI))
    print( dim( wPHI))}
# square root of precision matrix of the lattice process
#   solve(t(H)%*%H) is proportional to the covariance matrix of the Markov Random Field
  Q<-  LKrig.precision(LKinfo)

# variational matrix used to find coefficients of fitted surface and evaluate
  
############################################################################################
# this is the regularized regression matrix that is the key to the entire algorithm:
########################################################################################
  if(verbose){
    cat("spam class and dim for Q", fill=TRUE)
    print( is.spam(Q))
    print( dim( Q))}
  M<- t(wPHI)%*% wPHI + lambda*(Q)
  nzero <- length(M@entries)
#                            PC ERROR Message: Error in eval.with.vis(expr, envir, enclos) : 
#                           get slot "entries" from an object of a basic class ("matrix") with no slots
  if(verbose){
    cat("Number of nonzero elements:", nzero, fill=TRUE)}
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
#  the current "temp" matrix  has the same sparse pattern as that
#  which resulted in the factorization  cholesky as "use.cholesky"
  if(  is.null(use.cholesky)){
    Mc<- chol( M)}
  else{
    # reuse a previous decomposition to save computation
#    # first check that it might be from the same
#    if( sum( temp@colindices - use.cholesky@colindices)!=0){
#      stop("use.cholesky not the same sparse pattern as  t(wPHI)%*% wPHI + lambda*(Q) ")}
    Mc<-update.spam.chol.NgPeyton(use.cholesky, M)}
 
# If this is just to set up calculations return the intermediate results
# this list is used in MLE.LKrig
  out1<- LKrig.coef( Mc, wPHI, wT.matrix, wy, lambda)
  if( verbose){
    cat("d.coef", out1$d.coef, fill=TRUE)}
  fitted.values<- (wT.matrix%*%out1$d.coef + wPHI%*%out1$c.coef)/sqrt(weights)
  residuals<- y- fitted.values

  out2<-LKrig.lnPlike(Mc, Q, y,lambda,residuals, weights)
# save seed if random number generation happening outside LKrig
  if(exists(".Random.seed",1)){
    save.seed<- .Random.seed}
  else{
    save.seed <- NA} 
    if( !is.na(iseed)){
      set.seed(iseed)}
# generate N(0,1)  
  wEy<- matrix( rnorm( NtrA*n), n,NtrA)*sqrt(weights)
# restore seed   
   if( !is.na(iseed) & !is.na(save.seed[1])){
         assign(".Random.seed",save.seed, pos=1)}
#  
  out3<- LKrig.coef( Mc, wPHI, wT.matrix, wEy, lambda)
  wEyhat<-(wT.matrix%*%out3$d.coef + wPHI%*%out3$c.coef) 
  trA.info <- t((wEy *wEyhat)/weights)%*%rep(1,n)
  trA.est <- mean(trA.info)
  trA.SE  <- sd(trA.info)/ sqrt( length(trA.info))
  if(verbose){
    cat("trA.est", trA.est, fill=TRUE)}
# find the GCV function
  GCV=  (sum(weights*(residuals)^2)/n )  /( 1- trA.est/m)^2
   if(verbose){
    cat("GCV", GCV, fill=TRUE)}
# the output object
# note the ifelse switch whether to return the big cholesky decomposition  
  object<-list(x=x,y=y,weights=weights, Z=Z,
              d.coef=out1$d.coef, c.coef=out1$c.coef,
              fitted.values=fitted.values, residuals= residuals,
              LKinfo=LKinfo,
              GCV=GCV, lnProfileLike= out2$lnProfileLike,
              rho.MLE=out2$rho.MLE, shat.MLE= out2$shat.MLE,lambda=lambda,
              lnDetCov= out2$lnDetCov, quad.form= out2$quad.form,
              trA.info=trA.info, trA.est=trA.est, trA.SE=trA.SE,
              eff.df= trA.est, m=m, lambda.fixed=lambda,
              nonzero.entries=nzero,spatialdriftorder=2,nt=nt, 
              nZ=nZ,ind.drift=ind.drift,  # What is m? Now m_old?
              call=match.call())
  if( return.cholesky){
     object$Mc <-Mc}
# set the class and return.  
  class(object)<- "LKrig"
  return(object)
}

