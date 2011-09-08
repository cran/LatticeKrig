LKrig.coef <-
function(Mc,wPHI,wT.matrix,wy,lambda){
  if(length(lambda)>1){
    stop("lambda must be a scalar")}
  A <- forwardsolve(Mc, transpose = TRUE, t(wPHI)%*%wT.matrix, upper.tri = TRUE)
  A <- backsolve(Mc, A)
  A<- t(wT.matrix)%*%(wT.matrix - wPHI%*%A)/lambda
#   A is  (T^t M^{-1} T) 
  b <- forwardsolve(Mc, transpose = TRUE, t(wPHI)%*%wy, upper.tri = TRUE)
  b <- backsolve(Mc, b)  
  b<- t(wT.matrix)%*%(wy- wPHI%*%b)/lambda
# b is   (T^t M^{-1} y)
# Save the intermediate matrix   (T^t M^{-1} T) ^{-1}
# this the GLS covariance matrix of estimated coefficients
# should be small -- the default is 3X3  
  Omega<- solve( A)     # solve without right-hand argument calculates inverse
# GLS estimates   
  d.coef<- Omega%*% b
# coefficients of basis functions.   
  c.coef <- forwardsolve(Mc, transpose = TRUE, t(wPHI)%*%(wy-wT.matrix%*%d.coef),
                         upper.tri = TRUE)
  c.coef <- backsolve(Mc, c.coef)

  return(list( c.coef=c.coef, d.coef= d.coef))  
  }

