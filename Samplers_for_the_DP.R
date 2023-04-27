#####################################################
# Simulate from the DP - Polya Urn scheme
#####################################################
set.seed(54321)
# G_0 is a standard Normal
alphadp<-4 # the concentration parameter for the DP
N<-100 # number of samples
mu<-0 # normal mean for G_0
sigma<-1 # normal standard deviation for G_0
theta<-numeric(N)
indicator<-rep(0,N) # indicates the discrete part of the sampling with 1
theta[1]<-rnorm(1,mu,sigma)
for (i in 2:N){
  if (runif(1,0,1) < alphadp/(alphadp+(i-1))) {
    theta[i]<-rnorm(1,mu,sigma)} else {
      indicator[i]<-1
      forprobs<-1/(i-1)
      probs<-forprobs*c(rep(1,(i-1)))
      pick_theta<-rmultinom(1,1,probs) # vector with element 1 
      # where a theta is picked to be replicated
      indicator[which(pick_theta==1)]=1
      theta[i]<-sum(theta[1:(i-1)]*pick_theta)
    }
}

##########################################################################
# Simulate from the DP - Chinese restaurant process
##########################################################################
set.seed(54321)
# G_0 is a standard Normal
alphadp<-1000 # the concentration parameter for the DP
N<-100 # number of samples
mu<-0 # normal mean
sigma<-1 # normal standard deviation
theta<-numeric(N)
indicator<-rep(0,N) # indicates the discrete part of the sampling with 1
tableIndicator<-rep(0,N) # indicates the table an observation is assigned to
theta[1]<-rnorm(1,mu,sigma)
tableIndicator[1]<-1
for (i in 2:N){
  if (runif(1,0,1) < alphadp/(alphadp+(i-1))) {
    theta[i]<-rnorm(1,mu,sigma)
    tableIndicator[i]<-max(tableIndicator)+1
  } 
  else {
    indicator[i]<-1
    forprobs<-1/(i-1)
    probs<-forprobs*c(rep(1,(i-1))) 
    # vector of size (i-1) containing probabilities 1/(i-1)
    pick_theta<-rmultinom(1,1,probs) 
    # vector of size (i-1) with element 1 where a theta is picked to be replicated
    indicator[which(pick_theta==1)]=1
    theta[i]<-sum(theta[1:(i-1)]*pick_theta)
    tableIndicator[i]<-tableIndicator[which(pick_theta==1)]
  }
}

##########################################################################
# Simulate from the DP - Stick breaking representation
##########################################################################
set.seed(54321)
# G_0 is a standard Normal
alphadp<-1000 # the concentration parameter for the DP
N<-100 # number of samples
mu<-0 # normal mean
sigma<-1 # normal standard deviation
theta<-numeric(N)
indicator<-rep(0,N) 
# indicates the thetas that form the discrete part of the DP mean with 1
thetastarIndicator<-rep(0,N) 
# indicates the thetastar an observation is assigned to
v<-rep(1,N) # vector for the v quantities
theta[1]<-rnorm(1,mu,sigma)
thetastarIndicator[1]<-1
thetaAssigned<-0
v[1]<-rbeta(1,1,alphadp)
for (i in 2:N){
  maxthetastarIndicator<-max(thetastarIndicator)
  v[maxthetastarIndicator+1]<-rbeta(1,1,alphadp)
  psiProbs<-rep(0,maxthetastarIndicator+1)
  psiProbs[1]<-v[1]
  if (maxthetastarIndicator>1) {
    for (istar in 2:maxthetastarIndicator) {
      vForProduct<-c(1-v[1:(istar-1)],v[istar])
      psiProbs[istar]<-prod(vForProduct)
    }
  }
  psiProbs[maxthetastarIndicator+1]<-1-sum(psiProbs[1:maxthetastarIndicator])
  allocation<-which(as.vector(rmultinom(1,1,psiProbs)==1))
  thetastarIndicator[i]<-allocation
  if (thetastarIndicator[i]<(maxthetastarIndicator+1)) {
    theta[i]<-theta[which(thetastarIndicator==allocation)[1]]
    indicator[i]<-1
    indicator[which(thetastarIndicator==allocation)[1]]<-1
  } else {
    theta[i]<-rnorm(1,0,1)
    indicator[i]<-0
  }
}

