#####################################################################################
# TRUNCATED DPMM - EXAMPLE SAMPLER UNDER THE CRP REPRESENTATION (SECTION 4.1.1)
#####################################################################################
# FINITE K - RANDOM ALPHA - NEEDS AT LEAST 10000 BURNIN ITERATIONS TO CONVERGE 
#####################################################################################
set.seed(54321)
x<-matrix(c(0,0,0,
            0,0,0,
            0,0,0,
            0,0,0,
            0,0,0,
            1,1,1,
            1,1,1,
            1,1,1,
            1,1,1,
            1,1,1), nrow=10, byrow=TRUE)

iter<-100
burnin<-10000
rightLimitAlphadp<-6 # uniform prior for alphadp U(0.0001,rightLimitAlphadp)
alphag0<-0.5
betag0<-0.5
N<-10
K<-3
P<-3
z<-matrix(rep(0,iter*N), nrow=iter, ncol=N)
zcurrent<-matrix(rep(0,N), nrow=1, ncol=N)
znew<-matrix(rep(0,N), nrow=1, ncol=N)
phi <- array(0, dim = c(iter, K, P))
phicurrent <- array(0, dim = c(K, P))

# random group allocation to start 
for (isub in 1:N) {zcurrent[isub]<-which(as.vector(rmultinom(1,1,rep(1/K,K))==1))}
# phicurrent<-array(runif(K*P,0,1), dim = c(K, P)) # if we wanted to have a random \phi allocation to start
stDevForAlphadpProp<-0.2
alphadpcurrent<-1
alphadp<-rep(1,iter)

for (iiter in 2:(iter+burnin)) {
  #sample \phi
  for (icluster in 1:K) {
    for (ivar in 1:P) {
      numberOfSuccesses<-sum(x[zcurrent==icluster,ivar]==1)
      phicurrent[icluster,ivar]<-rbeta(1,alphag0+numberOfSuccesses,sum(zcurrent==icluster)+betag0-numberOfSuccesses)
    }
  }
  #sample alphadp
  alphadpproposed<-rnorm(1,alphadpcurrent,stDevForAlphadpProp)
  while ((alphadpproposed<0.0001) | (alphadpproposed>rightLimitAlphadp)) {
    alphadpproposed<-rnorm(1,alphadpcurrent,stDevForAlphadpProp)
  }
  MHrationumerator<-1/K
  for (iMHnum in 2:N) {
    MHrationumerator<-MHrationumerator*((sum(zcurrent[1:(iMHnum-1)]==zcurrent[iMHnum])+alphadpproposed/K)
                                        /((iMHnum-1)+alphadpproposed))
  }
  MHratiodenominator<-1/K
  for (iMHden in 2:N) {
    MHratiodenominator<-MHratiodenominator*((sum(zcurrent[1:(iMHden-1)]==zcurrent[iMHden])+alphadpcurrent/K)
                                            /((iMHden-1)+alphadpcurrent))
  }
  forTruncation1<-pnorm(rightLimitAlphadp,alphadpproposed,stDevForAlphadpProp)-pnorm(0.0001,alphadpproposed,stDevForAlphadpProp)
  forTruncation2<-pnorm(rightLimitAlphadp,alphadpcurrent,stDevForAlphadpProp)-pnorm(0.0001,alphadpcurrent,stDevForAlphadpProp)
  ratioOfQproposalsNum<-dnorm(alphadpcurrent,alphadpproposed,stDevForAlphadpProp)/forTruncation1
  ratioOfQproposalsDen<-dnorm(alphadpproposed,alphadpcurrent,stDevForAlphadpProp)/forTruncation2
  ratioOfQproposals<-ratioOfQproposalsNum/ratioOfQproposalsDen
  MHratio<-(MHrationumerator/MHratiodenominator) *ratioOfQproposals
  acceptProb<-min(c(1,MHratio))
  if (runif(1,0,1)<acceptProb) {alphadpcurrent<-alphadpproposed} 
  # sample z 
  for (isub in 1:N) {
    probForz<-rep(0,K)
    for (icluster in 1:K) {
      loglik<-sum(log(phicurrent[icluster,]^x[isub,])) + sum(log((1-phicurrent[icluster,])^(1-x[isub,])))
      nMinusic1<-0
      nMinusic2<-0
      if (isub!=1) {nMinusic1<-sum(znew[1:(isub-1)]==icluster)}
      if (isub!=N) {nMinusic2<-sum(zcurrent[(isub+1):N]==icluster)}
      nMinusic<-nMinusic1+nMinusic2
      probForz[icluster]<-exp(loglik)*(nMinusic+alphadpcurrent/K)/(N-1+alphadpcurrent)
    }
    probForz<-probForz/sum(probForz)
    znew[isub]<-which(as.vector(rmultinom(1,1,probForz)==1))
  }
  zcurrent<-znew
  # save samples
  if (iiter > burnin) {
    z[iiter-burnin,]<-znew
    phi[iiter-burnin,,]<-phicurrent
    alphadp[iiter-burnin]<-alphadpcurrent
  }
}
  
