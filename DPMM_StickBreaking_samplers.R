#################################################################################################
# DPMM - RANDOM ALPHA - WITHOUT INTEGRATING OUT \PSI - STICK BREAKING PRIOR FOR \PSI 
#################################################################################################
# UNIFORM PRIOR FOR ALPHA - GENERIC FULL CONDITIONAL FOR THE INFINITE MODEL (WALKER 2007)
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
burnin<-10000
rightLimitAlphadp<-8 # uniform prior for alphadp U(0.0001,rightLimitAlphadp)
lambdap<-0.5
N<-10
cMax<-14 # the maximum number of clusters just for building the arrays here
cInit<-4 # number of clusters for initial allocations
C<-10
P<-3
z<-matrix(rep(0,iter*N), nrow=iter, ncol=N)
zcurrent<-matrix(rep(0,N), nrow=1, ncol=N)
znew<-matrix(rep(0,N), nrow=1, ncol=N)
noOfLevels<-2
phi <- array(0, dim = c(iter, cMax, P,noOfLevels)) # for the StBreaking we sample the whole phi vector
phicurrent <- array(0, dim = c(cMax, P,noOfLevels))
psi<-array(0, dim = c(iter, cMax))
psicurrent<-matrix(rep(0,cMax), nrow=1, ncol=cMax)
vicurrent<-matrix(rep(0,cMax), nrow=1, ncol=cMax)
vicurrent[C]<-1
alphadp<-rep(1,iter)
alphadpcurrent<-1
ucurrent<-matrix(rep(0,N), nrow=1, ncol=N)
cStar<-cInit

# random group allocation to start 
for (isub in 1:N) {zcurrent[isub]<-which(as.vector(rmultinom(1,1,rep(1/cInit,cInit))==1))}
stDevForAlphadpProp<-0.02

for (iiter in 2:(iter+burnin)) {
  # sample \vi and calculate corresponding \psi
  clusterSizes<-rep(0,cStar)
  for (ics in 1:cStar) {clusterSizes[ics]<-sum(zcurrent==ics)}
  for (ivi in 1:(cStar-1)) {vicurrent[ivi]<-rbeta(1,clusterSizes[ivi]+1,sum(clusterSizes[(ivi+1):cStar])+alphadpcurrent)}
  psicurrent[1]<-vicurrent[1]
  for (ipsi in 2:cStar) {psicurrent[ipsi]<-vicurrent[ipsi]*prod(1-vicurrent[1:(ipsi-1)])}
  if ((sum(psicurrent<(10^(-10))))>0) {
    sumOfPsiZero<-sum(psicurrent<(10^(-10)))
    psicurrent[psicurrent<(10^(-10))]<-10^(-10)
    psicurrent[which(psicurrent==max(psicurrent))]<-psicurrent[which(psicurrent==max(psicurrent))]-sumOfPsiZero*(10^(-10))
  }
  #sample \phi
  for (icluster in 1:cStar) {
    for (ivar in 1:P) {
      numberOfSuccesses<-sum(x[zcurrent==icluster,ivar]==1)
      phicurrent[icluster,ivar,]<-rdirichlet(1,c(lambdap+numberOfSuccesses,sum(zcurrent==icluster)+lambdap-numberOfSuccesses))
    }
  }
  #sample alphadp
  alphadpproposed<-rnorm(1,alphadpcurrent,stDevForAlphadpProp)
  while ((alphadpproposed<0.0001) | (alphadpproposed>rightLimitAlphadp)) {
    alphadpproposed<-rnorm(1,alphadpcurrent,stDevForAlphadpProp)
  }
  dparam<-0 # the number of non-empty clusters
  for (iclust in 1:cStar) {if (sum(zcurrent==iclust)>0) {dparam<-dparam+1} }
  MHratio1<-(alphadpcurrent^dparam)*gamma(alphadpcurrent)/gamma(alphadpcurrent+N)
  MHratio2<-(alphadpproposed^dparam)*gamma(alphadpproposed)/gamma(alphadpproposed+N)
  forTruncation1<-pnorm(rightLimitAlphadp,alphadpproposed,stDevForAlphadpProp)-pnorm(0.0001,alphadpproposed,stDevForAlphadpProp)
  forTruncation2<-pnorm(rightLimitAlphadp,alphadpcurrent,stDevForAlphadpProp)-pnorm(0.0001,alphadpcurrent,stDevForAlphadpProp)
  ratioOfQproposalsNum<-dnorm(alphadpcurrent,alphadpproposed,stDevForAlphadpProp)/forTruncation1
  ratioOfQproposalsDen<-dnorm(alphadpproposed,alphadpcurrent,stDevForAlphadpProp)/forTruncation2
  ratioOfQproposals<-ratioOfQproposalsNum/ratioOfQproposalsDen
  MHratio<-(MHratio1/MHratio2) *ratioOfQproposals
  acceptProb<-min(c(1,MHratio))
  if (runif(1,0,1)<acceptProb) {alphadpcurrent<-alphadpproposed} 
  # sample u
  for (iu in 1:N) {ucurrent[iu]<-runif(1,0,psicurrent[zcurrent[iu]])}
  # find cstar, the required number of clusters
  cstar<-1
  while (sum(psicurrent[1:cstar])<=min(ucurrent)) {cstar<-cstar+1}
  # sample z 
  for (isub in 1:N) {
    probForz<-rep(0,cStar)
    for (icluster in 1:C) {
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

#################################################################################################
#################################################################################################
# DPMM - RANDOM ALPHA - NOT INTEGRATING OUT \PSI - STICK BREAKING PRIOR FOR \PSI 
#################################################################################################
# GAMMA PRIOR FOR ALPHA - USING THE FULL CONDITIONAL FOR THE INFINITE MODEL (ESCOBAR + WEST 1995)
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
burnin<-10000
rightLimitAlphadp<-8 # uniform prior for alphadp U(0.0001,rightLimitAlphadp)
gamma1<-0.5 # prior for alpha
gamma2<-0.5 # prior for alpha
lambdap<-0.5
N<-10
cMax<-14 # the maximum number of clusters just for building the arrays here
cInit<-4 # number of clusters for initial allocations
C<-10
P<-3
z<-matrix(rep(0,iter*N), nrow=iter, ncol=N)
zcurrent<-matrix(rep(0,N), nrow=1, ncol=N)
znew<-matrix(rep(0,N), nrow=1, ncol=N)
noOfLevels<-2
phi <- array(0, dim = c(iter, cMax, P,noOfLevels)) # for the StBreaking we sample the whole phi vector
phicurrent <- array(0, dim = c(cMax, P,noOfLevels))
psi<-array(0, dim = c(iter, cMax))
psicurrent<-matrix(rep(0,cMax), nrow=1, ncol=cMax)
vicurrent<-matrix(rep(0,cMax), nrow=1, ncol=cMax)
vicurrent[C]<-1
alphadp<-rep(1,iter)
alphadpcurrent<-1
ucurrent<-matrix(rep(0,N), nrow=1, ncol=N)
cStar<-cInit

# random group allocation to start 
for (isub in 1:N) {zcurrent[isub]<-which(as.vector(rmultinom(1,1,rep(1/cInit,cInit))==1))}
stDevForAlphadpProp<-0.02

for (iiter in 2:(iter+burnin)) {
  # sample \vi and calculate corresponding \psi
  clusterSizes<-rep(0,cStar)
  for (ics in 1:cStar) {clusterSizes[ics]<-sum(zcurrent==ics)}
  for (ivi in 1:(cStar-1)) {vicurrent[ivi]<-rbeta(1,clusterSizes[ivi]+1,sum(clusterSizes[(ivi+1):cStar])+alphadpcurrent)}
  psicurrent[1]<-vicurrent[1]
  for (ipsi in 2:cStar) {psicurrent[ipsi]<-vicurrent[ipsi]*prod(1-vicurrent[1:(ipsi-1)])}
  if ((sum(psicurrent<(10^(-10))))>0) {
    sumOfPsiZero<-sum(psicurrent<(10^(-10)))
    psicurrent[psicurrent<(10^(-10))]<-10^(-10)
    psicurrent[which(psicurrent==max(psicurrent))]<-psicurrent[which(psicurrent==max(psicurrent))]-sumOfPsiZero*(10^(-10))
  }
  #sample \phi
  for (icluster in 1:cStar) {
    for (ivar in 1:P) {
      numberOfSuccesses<-sum(x[zcurrent==icluster,ivar]==1)
      phicurrent[icluster,ivar,]<-rdirichlet(1,c(lambdap+numberOfSuccesses,sum(zcurrent==icluster)+lambdap-numberOfSuccesses))
    }
  }
  #sample alphadp
  dparam<-0 # the number of non-empty clusters
  for (iclust in 1:cStar) {if (sum(zcurrent==iclust)>0) {dparam<-dparam+1} }
  etaParam<-rbeta(1,alphadpcurrent+1,N)
  if (runif(1,0,1)<((gamma1+dparam-1)/((gamma1+dparam-1)+N*(gamma2-log(etaParam))))) {
    alphadpcurrent<-rgamma(1,gamma1+dparam,gamma2-log(etaParam))
  } 
  else {alphadpcurrent<-rgamma(1,gamma1+dparam-1,gamma2-log(etaParam))}
  # sample u
  for (iu in 1:N) {ucurrent[iu]<-runif(1,0,psicurrent[zcurrent[iu]])}
  # find cstar, the required number of clusters
  cstar<-1
  while (sum(psicurrent[1:cstar])<=min(ucurrent)) {cstar<-cstar+1}
  # sample z 
  for (isub in 1:N) {
    probForz<-rep(0,cStar)
    for (icluster in 1:C) {
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

