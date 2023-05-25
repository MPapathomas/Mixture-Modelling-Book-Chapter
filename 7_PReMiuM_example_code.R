########################################################################################################
# Continuous variables for the DPMM using PReMiuM (Section 5)
########################################################################################################
library(PReMiuM)
library(mvtnorm)
# Simulate data
set.seed(12345)
x <- matrix(0,500,4)
cov_matrix <- matrix(c(6,1.0,1.0,
                  1.0,6,1.0,
                  1.0,1.0,6), nrow=3, byrow=TRUE)
x[1:250,2:4] <- rmvnorm(250, mean=rep(32,3), sigma=cov_matrix)
x[251:500,2:4] <- rmvnorm(250, mean=rep(48,3), sigma=cov_matrix)
# Prepare data
xnames <- c('outcome', 'NO2', 'PM2_5', 'PM10')
x <- as.data.frame(x)
data_use <- x
colnames(data_use) <- paste0(xnames)
data_names <- c('NO2', 'PM2_5', 'PM10')
# Run the MCMC sampler
runInfoObj<-profRegr(excludeY=TRUE,  outcome="outcome", xModel="Normal", nBurn=5000, nSweeps=5000, 
                     seed=54321, nClusInit=30, data=data_use, covNames = data_names, 
                     output=paste("output",sep=""), nProgress=100, useNormInvWishPrior=TRUE, useHyperpriorR1=TRUE)
# Post-processing
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
clusObj$nClusters
clusObj$clusterSizes
clusObj$clustering
riskProfileObj<-calcAvgRiskAndProfile(clusObj, includeFixedEffects=F)
# Create summary figure
clusterOrderObj<-plotRiskProfile(riskProfileObj, "summary.png", showRelativeRisk=F, orderBy=NULL)

########################################################################################################
# Variable selection and the DPMM using PReMiuM (Section 6)
########################################################################################################
library(PReMiuM)
# Simulate data
set.seed(12345)
x <- matrix(0,500,6)
x[1:250,2:5] <- rbinom(1000,1,0.9)
x[251:500,2:5] <- rbinom(1000,1,0.1)
x[,6] <- rbinom(500,1,0.5)
# Prepare data
xnames <- c('LungCancer', 'NO2', 'PM2_5', 'PM10', 'Traffic', 'PhyAct')
x <- as.data.frame(x)
data_use <- x
colnames(data_use) <- paste0(xnames)
data_names <- c('NO2', 'PM2_5', 'PM10', 'Traffic', 'PhyAct')
# Run the MCMC sampler
runInfoObj<-profRegr(excludeY=TRUE,  outcome="outcome", xModel="Discrete", nBurn=5000, nSweeps=5000, 
                     seed=54321,  nClusInit=30, data=data_use, covNames = data_names, varSelectType ="BinaryCluster",
                     output=paste("output",sep=""), nProgress=100)
# Post-processing
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
clusObj$nClusters
clusObj$clusterSizes
clusObj$clustering
riskProfileObj<-calcAvgRiskAndProfile(clusObj, includeFixedEffects=F)
rho<-summariseVarSelectRho(runInfoObj)
# Choose variables with median selection probability > 0.5 for the summary figure 
# (All variables chosen in this toy example)
whichcovariates<-which(rho$rhoMedian>-0.1)
# Create summary figure
clusterOrderObj<-plotRiskProfile(riskProfileObj, "summary.png", whichCovariates=whichcovariates, showRelativeRisk=F, orderBy=NULL)

######################################################################################################
# Profile regression example code using PReMiuM (Section 7)
########################################################################################################
library(PReMiuM)
# Simulate data
set.seed(12345)
x <- matrix(0,500,6)
x[1:250,2:5] <- rbinom(1000,1,0.9)
x[251:500,2:5] <- rbinom(1000,1,0.1)
x[1:250,1] <- rbinom(250,1,0.7)
x[251:500,1] <- rbinom(250,1,0.1)
x[,6] <- rbinom(500,1,0.5)
# Prepare the data
xnames <- c('ShortnessBreath', 'NO2', 'PM2_5', 'PM10', 'Traffic', 'PhyAct')
x <- as.data.frame(x)
data_use <- x
colnames(data_use) <- paste0(xnames)
data_names <- c('NO2', 'PM2_5', 'PM10', 'Traffic', 'PhyAct')
# Run the MCMC sampler
runInfoObj<-profRegr(outcome=c("ShortnessBreath"),yModel="Bernoulli", xModel="Discrete", nSweeps=5000, 
                     nClusInit=30, seed=54321, varSelectType ="BinaryCluster",
                     nBurn=5000, data=data_use, output=paste("output",sep=""), nProgress=100, 
                     covNames = data_names)
# Post-processing
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
clusObj$nClusters
clusObj$clusterSizes
clusObj$clustering
riskProfileObj<-calcAvgRiskAndProfile(clusObj, includeFixedEffects=F)
rho<-summariseVarSelectRho(runInfoObj)
# Choose variables with median selection probability > 0.5 for the summary figure 
# (All variables chosen in this toy example)
whichcovariates<-which(rho$rhoMedian>-0.5)
# Change cluster labels according to the ordering of the clusters in terms of empirical risk
order_risk_ind <- order(riskProfileObj$empiricals)
# Choose clusters to be included in the summary figure (all clusters in this toy example)
whichclusters <- 1:clusObj$nClusters
# Create summary figure
clusterOrderObj<-plotRiskProfile(riskProfileObj, "summary.png", showRelativeRisk=F,orderBy="Empirical",
                                 whichCovariates=whichcovariates, whichClusters=whichclusters, useProfileStar=T)
