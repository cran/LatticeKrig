LKrig.lnPlike <-
function(Mc,Q,y,lambda,residuals,weights){
  n<- length(y)
  m<- dim(Q)[1]
  c.mKrig<- weights*residuals/lambda
# find log determinant of reg matrix temp for use in the log likeihood
  lnDetReg <- 2 * sum(log(diag(Mc)))
  lnDetQ<-  2* sum( log( diag( chol(Q))))
# now apply a miraculous determinant identity (Sylvester''s theorem)
#  det( I + UV) = det( I + VU)    providing UV is square
# or more generally   
#  det( A + UV) = det(A) det( I + V A^{-1}U)
#  or as we use it
#  ln(det( I + V A^{-1}U)) = ln( det(  A + UV)) - ln( det(A))
#  
#  with A = lambda* t(H)%*%H , U= t(PHI) V= PHI   
#  M =  lambda*(I + V A^{-1}U)
#
#   lnDet (M) = N* ln(lambda) + ln(Det(I + V A^{-1}U))
#             = N* ln(lambda) + ln(Det( temp)) - ln Det( lambda*Q)
#             = N* ln(lambda) + ln(Det( temp))  - ln Det( Q) - Ntotal* log( lambda)
# Here we  canceled out lambda terms 
# 
  lnDetCov<- lnDetReg - lnDetQ + (n - m)* log(lambda)
# finding quadratic form
# this uses a shortcut formula for the quadratic form in terms of
# the residuals 
  quad.form<-   sum(y* c.mKrig )
# MLE estimate of rho and sigma
# these are derived by assuming Y is  MN(  Td, rho*M )  
  rho.MLE <- quad.form/n
  shat.MLE <- sigma.MLE <- sqrt(lambda * rho.MLE)
# the  log profile likehood with  rhohat  and  dhat substituted
# leaving a profile for just lambda.
# note that this is _not_  -2*loglike just the log and
# includes the constants
  lnProfileLike <- (-n/2 - log(2*pi)*(n/2)
                      - (n/2)*log(rho.MLE) - (1/2) * lnDetCov)
  return( list(lnProfileLike=lnProfileLike,rho.MLE=rho.MLE, shat.MLE=shat.MLE,
               quad.form=quad.form, lnDetCov=lnDetCov) )
 
}

