library(PReMiuM)
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
dat_name <- c('NO2', 'PM2_5', 'PM10')
data_use
runInfoObj<-profRegr(excludeY=TRUE,  outcome="outcome", xModel="Discrete", 
                     nSweeps=5000, nClusInit=30, seed=54321, varSelectType ="BinaryCluster",
                     nBurn=5000, data=data_use, output=paste("output",sep=""), nProgress=100,
                     covNames = dat_name)
#post-processing
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
clusObj$nClusters
clusObj$clusterSizes
clusObj$clustering
riskProfileObj<-calcAvgRiskAndProfile(clusObj, includeFixedEffects=F)
rho<-summariseVarSelectRho(runInfoObj)
whichcovariates<-which(rho$rhoMedian>0.5)
whichcovariates
clusterOrderObj<-plotRiskProfile(riskProfileObj, "summary.png", showRelativeRisk=F, orderBy=NULL)

######################################################################################################

library(PReMiuM)
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
dat_name <- c('NO2', 'PM2_5', 'PM10')
data_use
runInfoObj<-profRegr(excludeY=TRUE,  outcome="outcome", xModel="Normal", nSweeps=5000, nClusInit=30, seed=54321,
                     nBurn=5000, data=data_use, output=paste("output",sep=""), nProgress=100,
                     covNames = dat_name)
#post-processing
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
clusObj$nClusters
clusObj$clusterSizes
clusObj$clustering
riskProfileObj<-calcAvgRiskAndProfile(clusObj, includeFixedEffects=F)
clusterOrderObj<-plotRiskProfile(riskProfileObj, "summary.png", showRelativeRisk=F, orderBy=NULL)

########################################################################################################

library(PReMiuM)
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
dat_name <- c('NO2', 'PM2_5', 'PM10')

runInfoObj<-profRegr(outcome=c("LungCancer"),yModel="Bernoulli", xModel="Discrete", nSweeps=5000, 
                     nClusInit=30, seed=54321, varSelectType ="BinaryCluster",
                     nBurn=5000, data=data_use, output=paste("output",sep=""), nProgress=100, 
                     covNames = dat_name)
#post-processing
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
clusObj$nClusters
clusObj$clusterSizes
clusObj$clustering
riskProfileObj<-calcAvgRiskAndProfile(clusObj, includeFixedEffects=F)
rho<-summariseVarSelectRho(runInfoObj)
whichcovariates<-which(rho$rhoMedian>0.5)
whichcovariates

dat_name[whichcovariates]
order_risk_ind <- order(riskProfileObj$empiricals) # cluster labels when ordering the clusters in terms of empirical risk
whichclusters <- 1:clusObj$nClusters
# need at most 4 covariates for plots to make sense
clusterOrderObj<-plotRiskProfile(riskProfileObj, "summary.png", showRelativeRisk=F,orderBy="Empirical",
                                 whichCovariates=whichcovariates, whichClusters=whichclusters, useProfileStar=T)
