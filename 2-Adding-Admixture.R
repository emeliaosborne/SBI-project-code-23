#Adding an Admixture Event.
#This R script contains code to perform ABC for model with 3 populations, 1 admixture event.

library("abc")
library("MASS")

sample_size<-20000
sample_meanvector <- c(0,0,0)  
mySigmaSq=1.25

#Function which generates a mvn sample and then summary statistics.
dataSimulation <- function(AmountOfData, sigma) {
  thisData<-mvrnorm(AmountOfData,
                    c(0,0,0), 
                    sigma)
  
  #Generate the summary statistic.
  thisCov<-cov(thisData)
  thisCovVec<-t(thisCov[upper.tri(thisCov,diag = T)])
  
  mySS<-data.frame(thisCovVec)
  
  return(mySS)
}

#Make the reference table.

simsSize<-20000

SSsimsX1Var = numeric(simsSize)
SSsimsCovX1X2 = numeric(simsSize)
SSsimsX2Var= numeric(simsSize)
SSsimsCovX1X3 = numeric(simsSize)
SSsimsCovX2X3 = numeric(simsSize)
SSsimsX3Var= numeric(simsSize)
SSsims<-data.frame(SSsimsX1Var,SSsimsCovX1X2,SSsimsX2Var,SSsimsCovX1X3, SSsimsCovX2X3, SSsimsX3Var)

delta1Sims<-rexp(simsSize,1/500)
deltamSims<-rexp(simsSize,1/500)
delta2MinusmSims<-rexp(simsSize,1/500)
alphaSims<-rbeta(simsSize,2,3)

paramSims<-data.frame(delta1Sims,deltamSims,delta2MinusmSims,alphaSims)
names(paramSims)<-c("delta1Sims", "deltamSims","delta2MinusmSims", "alphaSims")

for (val in 1:simsSize) {
  currentDelta1=delta1Sims[val]
  currentDeltam=deltamSims[val]
  currentDelta2Minusm=delta2MinusmSims[val]
  currentAlpha=alphaSims[val]
  currentCovMatrix=matrix(c(mySigmaSq*(currentDelta1+currentDeltam+currentDelta2Minusm),
                            currentAlpha*mySigmaSq*(currentDelta1+currentDeltam), 
                            0, 
                            currentAlpha*mySigmaSq*(currentDelta1+currentDeltam), 
                            (currentAlpha^2*(currentDelta1+currentDeltam)+(1-currentAlpha)^2*(currentDelta1+currentDeltam)+currentDelta2Minusm)*mySigmaSq, 
                            (1-currentAlpha)*mySigmaSq*currentDelta1, 
                            0, 
                            (1-currentAlpha)*mySigmaSq*currentDelta1, 
                            mySigmaSq*(currentDelta1+currentDeltam+currentDelta2Minusm)),nrow=3,ncol=3,byrow = TRUE)
  SSsims[val,]<-dataSimulation(sample_size,currentCovMatrix)
}

#Generate 'observed' data.

myDelta1=350
myDeltam=160
myDelta2Minusm=640
myAlpha=0.25

sample_covariance_matrix <- matrix(c(mySigmaSq*(myDelta1+myDeltam+myDelta2Minusm),
                                     myAlpha*mySigmaSq*(myDelta1+myDeltam), 
                                     0, 
                                     myAlpha*mySigmaSq*(myDelta1+myDeltam), 
                                     (myAlpha^2*(myDelta1+myDeltam)+(1-myAlpha)^2*(myDelta1+myDeltam)+myDelta2Minusm)*mySigmaSq, 
                                     (1-myAlpha)*mySigmaSq*myDelta1, 
                                     0, 
                                     (1-myAlpha)*mySigmaSq*myDelta1, 
                                     mySigmaSq*(myDelta1+myDeltam+myDelta2Minusm)),
                                   nrow=3,ncol=3,byrow = TRUE)
# create multivariate normal sample
sample_distribution <- mvrnorm(n = sample_size,
                               mu = sample_meanvector, 
                               Sigma = sample_covariance_matrix)

realDataCov<-cov(sample_distribution)
covs<-t(realDataCov[upper.tri(realDataCov,diag = T)])
dataSSs<-data.frame(covs)
names(dataSSs)<-c("var_X1","cov_X1_X2","var_X2","cov_X1_X3","cov_X2_X3","var_X3")

#Run the abc, currently 1% kept. 
#Rejection.
abcOutput = abc(target=dataSSs, param=paramSims, sumstat=SSsims,tol=0.01,method="rejection")
#Local linear regression adjustment.
abcOutputLLin = abc(target=dataSSs, param=paramSims, sumstat=SSsims,tol=0.01,method="loclinear")
#Neural network regression adjustment.
abcOutputNNet = abc(target=dataSSs, param=paramSims, sumstat=SSsims,tol=0.01,method="neuralnet")

