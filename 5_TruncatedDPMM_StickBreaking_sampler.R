#################################################################################################
# TRUNCATED DPMM  (SECTION 4.2)
#################################################################################################
# RANDOM ALPHADP - NOT INTEGRATING OUT \PSI - STICK BREAKING PRIOR FOR \PSI 
################################################################################################
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
burnin<-10000
rightLimitAlphadp<-3 # uniform prior for alphadp U(0.0001,rightLimitAlphadp)
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
vicurrent<-matrix(rep(0,K), nrow=1, ncol=K)
vicurrent[K]<-1
alphadp<-rep(1,iter)
alphadpcurrent<-1

# random group allocation to start 
for (isub in 1:N) {zcurrent[isub]<-which(as.vector(rmultinom(1,1,rep(1/K,K))==1))}
stDevForAlphadpProp<-0.02

for (iiter in 2:(iter+burnin)) {
  #sample \phi
  for (icluster in 1:K) {
    for (ivar in 1:P) {
      numberOfSuccesses<-sum(x[zcurrent==icluster,ivar]==1)
      phicurrent[icluster,ivar,]<-rdirichlet(1,c(lambdap+numberOfSuccesses,sum(zcurrent==icluster)+lambdap-numberOfSuccesses))
    }
  }
  # sample \vi and calculate corresponding \psi
  clusterSizes<-rep(0,K)
  for (ics in 1:K) {clusterSizes[ics]<-sum(zcurrent==ics)}
  for (ivi in 1:(K-1)) {vicurrent[ivi]<-rbeta(1,clusterSizes[ivi]+1,sum(clusterSizes[(ivi+1):K])+alphadpcurrent)}
  psicurrent[1]<-vicurrent[1]
  for (ipsi in 2:K) {psicurrent[ipsi]<-vicurrent[ipsi]*prod(1-vicurrent[1:(ipsi-1)])}
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
  MHratio1<-(beta(1,alphadpcurrent)/beta(1,alphadpproposed))^(K-1)
  MHratio2<-(prod(1-vicurrent[1:(K-1)]))^(alphadpproposed-alphadpcurrent)
  forTruncation1<-pnorm(rightLimitAlphadp,alphadpproposed,stDevForAlphadpProp)-pnorm(0.0001,alphadpproposed,stDevForAlphadpProp)
  forTruncation2<-pnorm(rightLimitAlphadp,alphadpcurrent,stDevForAlphadpProp)-pnorm(0.0001,alphadpcurrent,stDevForAlphadpProp)
  ratioOfQproposalsNum<-dnorm(alphadpcurrent,alphadpproposed,stDevForAlphadpProp)/forTruncation1
  ratioOfQproposalsDen<-dnorm(alphadpproposed,alphadpcurrent,stDevForAlphadpProp)/forTruncation2
  ratioOfQproposals<-ratioOfQproposalsNum/ratioOfQproposalsDen
  MHratio<-MHratio1*MHratio2 *ratioOfQproposals
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
