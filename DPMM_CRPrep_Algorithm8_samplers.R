##############################################################################################
# DPMM CRP - INFINITE K - ALGORITHM 8 FROM NEAL 2000 (SECTION 4.1.3)
##############################################################################################
# UNIFORM PRIOR ON ALPHA - USING THE ASYMPTOTIC EXPRESSION FOR THE TRUNCATED MODEL (NEAL 2000)
##############################################################################################
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
burnin<-20000
rightLimitAlphadp<-6 # uniform prior for alphadp U(0.0001,rightLimitAlphadp)
alphag0<-0.5
betag0<-0.5
N<-10
Cmax<-14 # the maximum number of clusters just for building the arrays here
Cinit<-4 # number of clusters for initial allocations
muAddedClusters<-3
P<-3
z<-matrix(rep(0,iter*N), nrow=iter, ncol=N)
zcurrent<-matrix(rep(0,N), nrow=1, ncol=N)
znew<-matrix(rep(0,N), nrow=1, ncol=N)
phi <- array(0, dim = c(iter, Cmax, P))
phicurrent <- array(0, dim = c(Cmax, P))

# random group allocation to start 
for (isub in 1:N) {zcurrent[isub]<-which(as.vector(rmultinom(1,1,rep(1/Cinit,Cinit))==1))}
# random \phi allocation to start
# phicurrent<-array(runif(Cmax*P,0,1), dim = c(Cmax, P)) # if we wanted random \phi allocation to start
stDevForAlphadpProp<-0.2
alphadpcurrent<-1
alphadp<-rep(1,iter)

for (iiter in 2:(iter+burnin)) {
  #########################################################################################
  #sample \phi 
  #################################################################################
  cBar<-max(zcurrent)
  for (icluster in 1:(cBar+muAddedClusters)) {
    for (ivar in 1:P) {
      numberOfSuccesses<-sum(x[zcurrent==icluster,ivar]==1)
      phicurrent[icluster,ivar]<-rbeta(1,alphag0+numberOfSuccesses,sum(zcurrent==icluster)+betag0-numberOfSuccesses)
    }
  }
  ##########################################################################################
  #sample alphadp
  ############################################################################################
  alphadpproposed<-rnorm(1,alphadpcurrent,stDevForAlphadpProp)
  while ((alphadpproposed<0.0001) | (alphadpproposed>rightLimitAlphadp)) {
    alphadpproposed<-rnorm(1,alphadpcurrent,stDevForAlphadpProp)
  }
  MHrationumerator<-1
  for (iMHnum in 2:N) {
    if (sum(zcurrent[1:(iMHnum-1)]==zcurrent[iMHnum])>0) {
      MHrationumerator<-MHrationumerator*(sum(zcurrent[1:(iMHnum-1)]==zcurrent[iMHnum])/((iMHnum-1)+alphadpproposed))
    }
  }
  MHratiodenominator<-1
  for (iMHden in 2:N) {
    if (sum(zcurrent[1:(iMHden-1)]==zcurrent[iMHden])>0) {
      MHratiodenominator<-MHratiodenominator*(sum(zcurrent[1:(iMHden-1)]==zcurrent[iMHden])/((iMHden-1)+alphadpcurrent))
    }
  }
  forTruncation1<-pnorm(rightLimitAlphadp,alphadpproposed,stDevForAlphadpProp)-pnorm(0.0001,alphadpproposed,stDevForAlphadpProp)
  forTruncation2<-pnorm(rightLimitAlphadp,alphadpcurrent,stDevForAlphadpProp)-pnorm(0.0001,alphadpcurrent,stDevForAlphadpProp)
  ratioOfQproposalsNum<-dnorm(alphadpcurrent,alphadpproposed,stDevForAlphadpProp)/forTruncation1
  ratioOfQproposalsDen<-dnorm(alphadpproposed,alphadpcurrent,stDevForAlphadpProp)/forTruncation2
  ratioOfQproposals<-ratioOfQproposalsNum/ratioOfQproposalsDen
  MHratio<-(MHrationumerator/MHratiodenominator) *ratioOfQproposals
  acceptProb<-min(c(1,MHratio))
  if (runif(1,0,1)<acceptProb) {alphadpcurrent<-alphadpproposed} 
  
  ################################################################################################# sample z 
  ##############################################################################################
  for (isub in 1:N) {
    # relabeling the z allocations 
    # relabel the whole vector z, including i, and then just check if z_i is on its own in a cluster and move it to cbar+1 with its phi
    zcurrentRelabel<-zcurrent 
    phicurrentRelabel<-phicurrent
    zlabel<-1
    for (iLabel in 1:max(zcurrent)) {
      if (sum(zcurrent==iLabel)>0) {
        zcurrentRelabel[zcurrent==iLabel]<-zlabel
        phicurrentRelabel[zlabel,]<-phicurrent[iLabel,]
        zlabel<-zlabel+1
      }
    }
    # if the subject i is on its own in a cluster and it is not at the maximum of the labeling, then put it at the maximum
    if ((sum(zcurrentRelabel==zcurrentRelabel[isub])==1) & (zcurrentRelabel[isub]<max(zcurrentRelabel))) { 
      maxztemp<-max(zcurrentRelabel)
      phitemp<-phicurrentRelabel[maxztemp,]
      zcurrentRelabel[zcurrentRelabel==maxztemp]<-zcurrentRelabel[isub]
      phicurrentRelabel[maxztemp,]<-phicurrentRelabel[zcurrentRelabel[isub],]
      zcurrentRelabel[isub]<-maxztemp
      phicurrentRelabel[zcurrentRelabel[isub],]<-phitemp
    }
    phicurrent<-phicurrentRelabel
    zcurrent<-zcurrentRelabel
    # if subject i is on its own in a cluster, then cBar is max(z)-1, otherwise it is max(z)
    if (sum(zcurrent==zcurrent[isub])==1)  {cBar<-max(zcurrent)-1} else {cBar<-max(zcurrent)}
    # sample the phi for the muAddedClusters empty additional clusters, when subject i is not alone in a cluster
    if (cBar==max(zcurrent)) {
      for (iPhicluster in (cBar+1):(cBar+muAddedClusters)) {
        phicurrent[iPhicluster,]<-rbeta(P,alphag0,betag0)
      }
    }
    # sample the phi for the muAddedClusters-1 empty additional clusters, when subject i is alone in a cluster
    if (cBar==(max(zcurrent)-1)) {
      for (iPhicluster in (cBar+1):(cBar+muAddedClusters)) {
        phicurrent[iPhicluster,]<-rbeta(P,alphag0,betag0)
      }
    }
    # calculate the probabilities for the different clusters
    probForz<-rep(0,cBar+muAddedClusters)
    for (icluster in 1:cBar) {
      loglik<-sum(log(phicurrent[icluster,]^x[isub,])) + sum(log((1-phicurrent[icluster,])^(1-x[isub,])))
      nMinusic1<-0
      nMinusic2<-0
      if (isub!=1) {nMinusic1<-sum(zcurrent[1:(isub-1)]==icluster)}
      if (isub!=N) {nMinusic2<-sum(zcurrent[(isub+1):N]==icluster)}
      nMinusic<-nMinusic1+nMinusic2
      probForz[icluster]<-exp(loglik)*(nMinusic)/(N-1+alphadpcurrent)
    }
    for (icluster in (cBar+1):(cBar+muAddedClusters)) {
      loglik<-sum(log(phicurrent[icluster,]^x[isub,])) + sum(log((1-phicurrent[icluster,])^(1-x[isub,])))
      nMinusic1<-0
      nMinusic2<-0
      if (isub!=1) {nMinusic1<-sum(zcurrent[1:(isub-1)]==icluster)}
      if (isub!=N) {nMinusic2<-sum(zcurrent[(isub+1):N]==icluster)}
      nMinusic<-nMinusic1+nMinusic2
      probForz[icluster]<-exp(loglik)*(alphadpcurrent/muAddedClusters)/(N-1+alphadpcurrent)
    }
    probForz<-probForz/sum(probForz)
    zcurrent[isub]<-which(as.vector(rmultinom(1,1,probForz)==1))
  }
  
  # save samples
  if (iiter > burnin) {
    z[iiter-burnin,]<-zcurrent
    phi[iiter-burnin,,]<-phicurrent
    alphadp[iiter-burnin]<-alphadpcurrent
  }
}

##############################################################################################
##############################################################################################
# DPMM CRP - INFINITE K - ALGORITHM 8 FROM NEAL 2000 (SECTION 4.1.3)
##############################################################################################
# GAMMA PRIOR FOR ALPHA - USING THE FULL CONDITIONAL FOR THE INFINITE MODEL (ESCOBAR + WEST 1995)
##############################################################################################
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
burnin<-20000
rightLimitAlphadp<-6 # uniform prior for alphadp U(0.0001,rightLimitAlphadp)
alphag0<-0.5
betag0<-0.5
gamma1<-0.5 # prior for alpha
gamma2<-0.5 # prior for alpha
N<-10
Cmax<-14 # the maximum number of clusters just for building the arrays here
Cinit<-4 # number of clusters for initial allocations
muAddedClusters<-3
P<-3
z<-matrix(rep(0,iter*N), nrow=iter, ncol=N)
zcurrent<-matrix(rep(0,N), nrow=1, ncol=N)
znew<-matrix(rep(0,N), nrow=1, ncol=N)
phi <- array(0, dim = c(iter, Cmax, P))
phicurrent <- array(0, dim = c(Cmax, P))

# random group allocation to start 
for (isub in 1:N) {zcurrent[isub]<-which(as.vector(rmultinom(1,1,rep(1/Cinit,Cinit))==1))}
# random \phi allocation to start
# phicurrent<-array(runif(Cmax*P,0,1), dim = c(Cmax, P)) # if we wanted random \phi allocation to start
stDevForAlphadpProp<-0.2
alphadpcurrent<-1
alphadp<-rep(1,iter)

for (iiter in 2:(iter+burnin)) {
  ########################################################################################
  #sample \phi
  #########################################################################################
  cBar<-max(zcurrent)
  for (icluster in 1:(cBar+muAddedClusters)) {
    for (ivar in 1:P) {
      numberOfSuccesses<-sum(x[zcurrent==icluster,ivar]==1)
      phicurrent[icluster,ivar]<-rbeta(1,alphag0+numberOfSuccesses,sum(zcurrent==icluster)+betag0-numberOfSuccesses)
    }
  }
  ###########################################################################################
  #sample alphadp
  #################################################################################################
  dparam<-0 # the number of non-empty clusters
  for (iclust in 1:Cmax) {if (sum(zcurrent==iclust)>0) {dparam<-dparam+1} }
  etaParam<-rbeta(1,alphadpcurrent+1,N)
  if (runif(1,0,1)<((gamma1+dparam-1)/((gamma1+dparam-1)+N*(gamma2-log(etaParam))))) {
    alphadpcurrent<-rgamma(1,gamma1+dparam,gamma2-log(etaParam))
  } 
  else {alphadpcurrent<-rgamma(1,gamma1+dparam-1,gamma2-log(etaParam))}
  
  ###########################################################################################
  # sample z 
  ###############################################################################################
  for (isub in 1:N) {
    # relabeling the z allocations 
    # relabel the whole vector z, including i, and then just check if z_i is on its own in a cluster and move it to cbar+1 with its phi
    zcurrentRelabel<-zcurrent 
    phicurrentRelabel<-phicurrent
    zlabel<-1
    for (iLabel in 1:max(zcurrent)) {
      if (sum(zcurrent==iLabel)>0) {
        zcurrentRelabel[zcurrent==iLabel]<-zlabel
        phicurrentRelabel[zlabel,]<-phicurrent[iLabel,]
        zlabel<-zlabel+1
      }
    }
    # if the subject i is on its own in a cluster and it is not at the maximum of the labeling, then put it at the maximum
    if ((sum(zcurrentRelabel==zcurrentRelabel[isub])==1) & (zcurrentRelabel[isub]<max(zcurrentRelabel))) { 
      maxztemp<-max(zcurrentRelabel)
      phitemp<-phicurrentRelabel[maxztemp,]
      zcurrentRelabel[zcurrentRelabel==maxztemp]<-zcurrentRelabel[isub]
      phicurrentRelabel[maxztemp,]<-phicurrentRelabel[zcurrentRelabel[isub],]
      zcurrentRelabel[isub]<-maxztemp
      phicurrentRelabel[zcurrentRelabel[isub],]<-phitemp
    }
    phicurrent<-phicurrentRelabel
    zcurrent<-zcurrentRelabel
    # if subject i is on its own in a cluster, then cBar is max(z)-1, otherwise it is max(z)
    if (sum(zcurrent==zcurrent[isub])==1)  {cBar<-max(zcurrent)-1} else {cBar<-max(zcurrent)}
    # sample the phi for the muAddedClusters empty additional clusters, when subject i is not alone in a cluster
    if (cBar==max(zcurrent)) {
      for (iPhicluster in (cBar+1):(cBar+muAddedClusters)) {
        phicurrent[iPhicluster,]<-rbeta(P,alphag0,betag0)
      }
    }
    # sample the phi for the muAddedClusters-1 empty additional clusters, when subject i is alone in a cluster
    if (cBar==(max(zcurrent)-1)) {
      for (iPhicluster in (cBar+1):(cBar+muAddedClusters)) {
        phicurrent[iPhicluster,]<-rbeta(P,alphag0,betag0)
      }
    }
    # calculate the probabilities for the different clusters
    probForz<-rep(0,cBar+muAddedClusters)
    for (icluster in 1:cBar) {
      loglik<-sum(log(phicurrent[icluster,]^x[isub,])) + sum(log((1-phicurrent[icluster,])^(1-x[isub,])))
      nMinusic1<-0
      nMinusic2<-0
      if (isub!=1) {nMinusic1<-sum(zcurrent[1:(isub-1)]==icluster)}
      if (isub!=N) {nMinusic2<-sum(zcurrent[(isub+1):N]==icluster)}
      nMinusic<-nMinusic1+nMinusic2
      probForz[icluster]<-exp(loglik)*(nMinusic)/(N-1+alphadpcurrent)
    }
    for (icluster in (cBar+1):(cBar+muAddedClusters)) {
      loglik<-sum(log(phicurrent[icluster,]^x[isub,])) + sum(log((1-phicurrent[icluster,])^(1-x[isub,])))
      nMinusic1<-0
      nMinusic2<-0
      if (isub!=1) {nMinusic1<-sum(zcurrent[1:(isub-1)]==icluster)}
      if (isub!=N) {nMinusic2<-sum(zcurrent[(isub+1):N]==icluster)}
      nMinusic<-nMinusic1+nMinusic2
      probForz[icluster]<-exp(loglik)*(alphadpcurrent/muAddedClusters)/(N-1+alphadpcurrent)
    }
    probForz<-probForz/sum(probForz)
    zcurrent[isub]<-which(as.vector(rmultinom(1,1,probForz)==1))
  }
  
  # save samples
  if (iiter > burnin) {
    z[iiter-burnin,]<-zcurrent
    phi[iiter-burnin,,]<-phicurrent
    alphadp[iiter-burnin]<-alphadpcurrent
  }
}


