#################################################################################################
# SAMPLER FOR THE FINITE MIXTURE MODEL (SECTION 3.2)
#################################################################################################
# FINITE K - RANDOM ALPHA - NOT INTEGRATING OUT \PSI - DIRICHLET PRIOR FOR \PSI 
#################################################################################################
library(gtools)
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
iter<-1000
burnin<-5000
rightLimitAlphadp<-6 # uniform prior for alphadp U(0.0001,rightLimitAlphadp)
lambdap<-0.5
N<-10
K<-10
P<-3
z<-matrix(rep(0,iter*N), nrow=iter, ncol=N)
zcurrent<-matrix(rep(0,N), nrow=1, ncol=N)
znew<-matrix(rep(0,N), nrow=1, ncol=N)
noOfLevels<-2
phi <- array(0, dim = c(iter, K, P,noOfLevels)) # for the StBreaking we sample the whole phi vector
phicurrent <- array(0, dim = c(K, P,noOfLevels))
psi<-array(0, dim = c(iter, K))
psicurrent<-matrix(rep(0,K), nrow=1, ncol=K)
alphadp<-rep(1,iter)
alphadpcurrent<-1

# random group allocation to start 
for (isub in 1:N) {zcurrent[isub]<-which(as.vector(rmultinom(1,1,rep(1/K,K))==1))}
stDevForAlphadpProp<-0.2

for (iiter in 2:(iter+burnin)) {
  print(iiter)
  #sample \phi
  for (icluster in 1:K) {
    for (ivar in 1:P) {
      numberOfSuccesses<-sum(x[zcurrent==icluster,ivar]==1)
      phicurrent[icluster,ivar,]<-rdirichlet(1,c(lambdap+numberOfSuccesses,sum(zcurrent==icluster)+lambdap-numberOfSuccesses))
    }
  }
  # sample \psi
  forPsiProbs<-sum(zcurrent==1)
  for (ipsiprobs in 2:K) {forPsiProbs<-c(forPsiProbs,sum(zcurrent==ipsiprobs))}
  psicurrent<-rdirichlet(1,forPsiProbs+alphadpcurrent)
  if ((sum(psicurrent<(10^(-10))))>0) {
    sumOfPsiZero<-sum(psicurrent<(10^(-10)))
    psicurrent[psicurrent<(10^(-10))]<-10^(-10)
    psicurrent[which(psicurrent==max(psicurrent))]<-psicurrent[which(psicurrent==max(psicurrent))]-sumOfPsiZero*(10^(-10))
  }
  #sample alphadp
  alphadpproposed<-rnorm(1,alphadpcurrent,stDevForAlphadpProp)
  while ((alphadpproposed<0.0001) | (alphadpproposed>rightLimitAlphadp)) {
    alphadpproposed<-rnorm(1,alphadpcurrent,stDevForAlphadpProp)
  }
  MHrationumerator<-sum(log(psicurrent^(alphadpproposed-1)))
  MHratiodenominator<-sum(log(psicurrent^(alphadpcurrent-1)))
  forTruncation1<-pnorm(rightLimitAlphadp,alphadpproposed,stDevForAlphadpProp)-pnorm(0.0001,alphadpproposed,stDevForAlphadpProp)
  forTruncation2<-pnorm(rightLimitAlphadp,alphadpcurrent,stDevForAlphadpProp)-pnorm(0.0001,alphadpcurrent,stDevForAlphadpProp)
  ratioOfQproposalsNum<-dnorm(alphadpcurrent,alphadpproposed,stDevForAlphadpProp)/forTruncation1
  ratioOfQproposalsDen<-dnorm(alphadpproposed,alphadpcurrent,stDevForAlphadpProp)/forTruncation2
  ratioOfQproposals<-ratioOfQproposalsNum/ratioOfQproposalsDen
  MHratio<-exp(MHrationumerator-MHratiodenominator) *ratioOfQproposals
  acceptProb<-min(c(1,MHratio))
  if (runif(1,0,1)<acceptProb) {alphadpcurrent<-alphadpproposed} 
  # sample z 
  for (isub in 1:N) {
    probForz<-rep(0,K)
    for (icluster in 1:K) {
      probForz[icluster]<-psicurrent[icluster] 
      for (iForProbForz in 1:P) {probForz[icluster]<-probForz[icluster]*phicurrent[icluster,iForProbForz,x[isub,iForProbForz]+1] }
    }
    probForz<-probForz/sum(probForz)
    znew[isub]<-which(as.vector(rmultinom(1,1,probForz)==1))
  }
  zcurrent<-znew
  # save samples
  if (iiter > burnin) {
    z[iiter-burnin,]<-zcurrent
    phi[iiter-burnin,,,]<-phicurrent
    psi[iiter-burnin,]<-psicurrent
    alphadp[iiter-burnin]<-alphadpcurrent
  }
}
