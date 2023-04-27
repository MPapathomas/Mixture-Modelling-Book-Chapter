########################################################################################################
# Continuous variables for the DPMM using PReMiuM (Section 5)
########################################################################################################
library(PReMiuM)
# Prepare data
xnames <- c('outcome', 'NO2', 'PM2_5', 'PM10')
x<-matrix(c(0,31,37,24,
            0,33,32,26,
            0,26,34,29,
            0,28,33,32,
            0,30,25,29,
            1,48,55,44,
            1,53,50,43,
            1,46,49,48,
            1,53,53,51,
            1,52,52,48), nrow=10, byrow=TRUE)
x <- as.data.frame(x)
data_use <- x
colnames(data_use) <- paste0(xnames)
data_names <- c('NO2', 'PM2_5', 'PM10')
# Run the MCMC sampler
runInfoObj<-profRegr(excludeY=TRUE,  outcome="outcome", xModel="Normal", nBurn=5000, nSweeps=5000, 
                     seed=54321, nClusInit=30, data=data_use, covNames = data_names, 
                     output=paste("output",sep=""), nProgress=100)
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
# Prepare data
xnames <- c('LungCancer', 'NO2', 'PM2_5', 'PM10')
x<-matrix(c(0,0,0,0,
            0,0,0,0,
            0,0,0,0,
            0,0,0,0,
            0,0,0,0,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1), nrow=10, byrow=TRUE)
x <- as.data.frame(x)
data_use <- x
colnames(data_use) <- paste0(xnames)
data_names <- c('NO2', 'PM2_5', 'PM10')
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
whichcovariates<-which(rho$rhoMedian>0.5)
# Create summary figure
clusterOrderObj<-plotRiskProfile(riskProfileObj, "summary.png", whichCovariates=whichcovariates, showRelativeRisk=F, orderBy=NULL)

######################################################################################################
# Profile regression example code using PReMiuM (Section 7)
########################################################################################################
library(PReMiuM)
# Prepare the data
xnames <- c('LungCancer', 'NO2', 'PM2_5', 'PM10')
x<-matrix(c(0,0,0,0,
            0,0,0,0,
            0,0,0,0,
            0,0,0,0,
            0,0,0,0,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1,
            1,1,1,1), nrow=10, byrow=TRUE)
x <- as.data.frame(x)
data_use <- x
colnames(data_use) <- paste0(xnames)
data_names <- c('NO2', 'PM2_5', 'PM10')
# Run the MCMC sampler
runInfoObj<-profRegr(outcome=c("LungCancer"),yModel="Bernoulli", xModel="Discrete", nSweeps=5000, 
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
whichcovariates<-which(rho$rhoMedian>0.5)
# Change cluster labels according to the ordering of the clusters in terms of empirical risk
order_risk_ind <- order(riskProfileObj$empiricals)
# Choose clusters to be included in the summary figure (all clusters in this toy example)
whichclusters <- 1:clusObj$nClusters
# Create summary figure
clusterOrderObj<-plotRiskProfile(riskProfileObj, "summary.png", showRelativeRisk=F,orderBy="Empirical",
                                 whichCovariates=whichcovariates, whichClusters=whichclusters, useProfileStar=T)
