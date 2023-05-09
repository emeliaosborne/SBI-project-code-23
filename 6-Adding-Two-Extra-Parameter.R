#Adding an 2 extra parameters to the topology.
#This R script contains code to perform ABC for model with 3 populations, 2 admixture events, 6 parameters.
#Uses two sources of data, relationship between known.

library("abc")
library("MASS")

sample_size<-10000
sample_meanvector <- c(0,0,0) 
mySigmaSq=1.25
myLambda=1/500

#Functions which generates mvn samples and then summary statistics (1 and 2 sources).

dataSimulationN <- function(AmountOfData, sigma) {
  thisData<-mvrnorm(AmountOfData,
                    c(0,0,0), 
                    sigma)
  
  #Generate the summary statistic.
  thisCov<-cov(thisData)
  thisCovVec<-t(thisCov[upper.tri(thisCov,diag = T)])
  
  mySS<-data.frame(thisCovVec)
  
  return(mySS)
}

dataSimulation <- function(AmountOfData, sigma, sigmaB) {
  thisData<-mvrnorm(AmountOfData,
                    c(0,0,0), 
                    sigma)
  
  thisDataB<-mvrnorm(AmountOfData,
                     c(0,0,0), 
                     sigmaB)
  
  #Generate the summary statistic.
  thisCov<-cov(thisData)
  thisCovB<-cov(thisDataB)
  thisCovVec<-t(thisCov[upper.tri(thisCov,diag = T)])
  thisCovVecB<-t(thisCovB[upper.tri(thisCovB,diag = T)])
  
  mySS<-data.frame(thisCovVec,thisCovVecB)
  
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

SSsimsX1VarB = numeric(simsSize)
SSsimsCovX1X2B = numeric(simsSize)
SSsimsX2VarB= numeric(simsSize)
SSsimsCovX1X3B = numeric(simsSize)
SSsimsCovX2X3B = numeric(simsSize)
SSsimsX3VarB= numeric(simsSize)
SSsims<-data.frame(SSsimsX1Var,SSsimsCovX1X2,SSsimsX2Var,SSsimsCovX1X3, SSsimsCovX2X3, SSsimsX3Var,SSsimsX1VarB,SSsimsCovX1X2B,SSsimsX2VarB,SSsimsCovX1X3B, SSsimsCovX2X3B, SSsimsX3VarB)

delta1Sims<-rexp(simsSize,1/500)
deltamSims<-rexp(simsSize,1/500)
deltacSims<-rexp(simsSize,1/500)
delta2MinusmMinuscSims<-rexp(simsSize,1/500)
alphaSims<-rbeta(simsSize,2,3)
betaSims<-rbeta(simsSize,2,3)

paramSims<-data.frame(delta1Sims,deltamSims,deltacSims,delta2MinusmMinuscSims,alphaSims,betaSims)
names(paramSims)<-c("delta1Sims", "deltamSims", "deltacSims","delta2MinusmMinuscSims", "alphaSims","betaSims")

for (val in 1:simsSize) {
  currentDelta1=delta1Sims[val]
  currentDeltam=deltamSims[val]
  currentDeltac=deltacSims[val]
  currentDelta2MinusmMinusc=delta2MinusmMinuscSims[val]
  currentAlpha=alphaSims[val]
  currentBeta=betaSims[val]
  
  Var1C=mySigmaSq*(currentAlpha^2*(currentDelta1+currentDeltam)+(1-currentAlpha)^2*(currentDelta1+currentDeltam)+currentDeltac)
  CovX1X2=currentAlpha*mySigmaSq*(currentDelta1+currentDeltam)
  CovX1X3=currentAlpha*currentBeta*mySigmaSq*(currentDelta1+currentDeltam)
  CovX2X3=(1-currentAlpha)*(1-currentBeta)*mySigmaSq*currentDelta1+currentBeta*Var1C
  
  VarX2B=mySigmaSq/(2*myLambda)*(
    currentAlpha^2*(exp(2*myLambda*(currentDelta1+currentDeltam))-1)+
      (1-currentAlpha)^2*(exp(2*myLambda*(currentDelta1+currentDeltam))-1)+
      exp(2*myLambda*(currentDelta1+currentDeltam+currentDeltac+currentDelta2MinusmMinusc))-
      exp(2*myLambda*(currentDelta1+currentDeltam)))
  
  Var1CB=mySigmaSq/(2*myLambda)*(
    currentAlpha^2*(exp(2*myLambda*(currentDelta1+currentDeltam))-1)+
      (1-currentAlpha)^2*(exp(2*myLambda*(currentDelta1+currentDeltam))-1)+
      exp(2*myLambda*(currentDelta1+currentDeltam+currentDeltac))-
      exp(2*myLambda*(currentDelta1+currentDeltam)))
  CovX1X2B=currentAlpha*mySigmaSq/(2*myLambda)*(exp(2*myLambda*(currentDelta1+currentDeltam))-1)
  CovX1X3B=currentAlpha*currentBeta*mySigmaSq/(2*myLambda)*(exp(2*myLambda*(currentDelta1+currentDeltam))-1)
  CovX2X3B=(1-currentAlpha)*(1-currentBeta)*mySigmaSq/(2*myLambda)*(exp(2*myLambda*currentDelta1)-1)+currentBeta*Var1CB
  
  
  currentCovMatrix <- matrix(c(mySigmaSq*(currentDelta1+currentDeltam+currentDeltac+currentDelta2MinusmMinusc),
                               CovX1X2, 
                               CovX1X3, 
                               CovX1X2, 
                               mySigmaSq*(currentAlpha^2*(currentDelta1+currentDeltam)+(1-currentAlpha)^2*(currentDelta1+currentDeltam)+currentDeltac+currentDelta2MinusmMinusc), 
                               CovX2X3, 
                               CovX1X3, 
                               CovX2X3, 
                               currentBeta^2*Var1C+mySigmaSq*((1-currentBeta)^2*(currentDelta1+currentDeltam+currentDeltac+currentDelta2MinusmMinusc)+2*currentBeta*(1-currentBeta)*(1-currentAlpha)*currentDelta1+currentDelta2MinusmMinusc)),nrow=3,ncol=3,byrow = TRUE)
  
  currentCovMatrixB <- matrix(c(mySigmaSq/(2*myLambda)*(exp(2*myLambda*(currentDelta1+currentDeltam+currentDeltac+currentDelta2MinusmMinusc))-1),
                                CovX1X2B, 
                                CovX1X3B, 
                                CovX1X2B, 
                                VarX2B, 
                                CovX2X3B, 
                                CovX1X3B, 
                                CovX2X3B, 
                                currentBeta^2*Var1CB+mySigmaSq/(2*myLambda)*((1-currentBeta)^2*(exp(2*myLambda*(currentDelta1+currentDeltam+currentDeltac))-1)+2*currentBeta*(1-currentBeta)*(1-currentAlpha)*(exp(2*myLambda*currentDelta1)-1)+exp(2*myLambda*(currentDelta1+currentDeltam+currentDeltac+currentDelta2MinusmMinusc))-exp(2*myLambda*(currentDelta1+currentDeltam+currentDeltac)))
  ),nrow=3,ncol=3,byrow = TRUE)
  
  SSsims[val,]<-dataSimulation(sample_size,currentCovMatrix,currentCovMatrixB)
}
#Generate 'observed' data.

myDelta1=350
myDeltam=160
myDeltac=300
myDelta2MinusmMinusc=340
myAlpha=0.25
mySigmaSq=1.25
myBeta=0.4

Var1C=mySigmaSq*(myAlpha^2*(myDelta1+myDeltam)+(1-myAlpha)^2*(myDelta1+myDeltam)+myDeltac)
CovX1X2=myAlpha*mySigmaSq*(myDelta1+myDeltam)
CovX1X3=myAlpha*myBeta*mySigmaSq*(myDelta1+myDeltam)
CovX2X3=(1-myAlpha)*(1-myBeta)*mySigmaSq*myDelta1+myBeta*Var1C


VarX2B=mySigmaSq/(2*myLambda)*(
  myAlpha^2*(exp(2*myLambda*(myDelta1+myDeltam))-1)+
    (1-myAlpha)^2*(exp(2*myLambda*(myDelta1+myDeltam))-1)+
    exp(2*myLambda*(myDelta1+myDeltam+myDeltac+myDelta2MinusmMinusc))-
    exp(2*myLambda*(myDelta1+myDeltam)))

Var1CB=mySigmaSq/(2*myLambda)*(
  myAlpha^2*(exp(2*myLambda*(myDelta1+myDeltam))-1)+
    (1-myAlpha)^2*(exp(2*myLambda*(myDelta1+myDeltam))-1)+
    exp(2*myLambda*(myDelta1+myDeltam+myDeltac))-
    exp(2*myLambda*(myDelta1+myDeltam)))
CovX1X2B=myAlpha*mySigmaSq/(2*myLambda)*(exp(2*myLambda*(myDelta1+myDeltam))-1)
CovX1X3B=myAlpha*myBeta*mySigmaSq/(2*myLambda)*(exp(2*myLambda*(myDelta1+myDeltam))-1)
CovX2X3B=(1-myAlpha)*(1-myBeta)*mySigmaSq/(2*myLambda)*(exp(2*myLambda*myDelta1)-1)+myBeta*Var1CB


sample_covariance_matrix <- matrix(c(mySigmaSq*(myDelta1+myDeltam+myDeltac+myDelta2MinusmMinusc),
                                     CovX1X2, 
                                     CovX1X3, 
                                     CovX1X2, 
                                     mySigmaSq*(myAlpha^2*(myDelta1+myDeltam)+(1-myAlpha)^2*(myDelta1+myDeltam)+myDeltac+myDelta2MinusmMinusc), 
                                     CovX2X3, 
                                     CovX1X3, 
                                     CovX2X3, 
                                     myBeta^2*Var1C+mySigmaSq*((1-myBeta)^2*(myDelta1+myDeltam+myDeltac+myDelta2MinusmMinusc)+2*myBeta*(1-myBeta)*(1-myAlpha)*myDelta1+myDelta2MinusmMinusc)),nrow=3,ncol=3,byrow = TRUE)

sample_covariance_matrixB <- matrix(c(mySigmaSq/(2*myLambda)*(exp(2*myLambda*(myDelta1+myDeltam+myDeltac+myDelta2MinusmMinusc))-1),
                                      CovX1X2B, 
                                      CovX1X3B, 
                                      CovX1X2B, 
                                      VarX2B, 
                                      CovX2X3B, 
                                      CovX1X3B, 
                                      CovX2X3B, 
                                      myBeta^2*Var1CB + mySigmaSq/(2*myLambda)*((1-myBeta)^2*(exp(2*myLambda*(myDelta1+myDeltam+myDeltac))-1)+2*myBeta*(1-myBeta)*(1-myAlpha)*(exp(2*myLambda*myDelta1)-1)+exp(2*myLambda*(myDelta1+myDeltam+myDeltac+myDelta2MinusmMinusc))-exp(2*myLambda*(myDelta1+myDeltam+myDeltac)))),nrow=3,ncol=3,byrow = TRUE)

# create 2 multivariate normal samples
sample_distribution <- mvrnorm(n = sample_size,
                               mu = sample_meanvector, 
                               Sigma = sample_covariance_matrix)

sample_distributionB <- mvrnorm(n = sample_size,
                                mu = sample_meanvector, 
                                Sigma = sample_covariance_matrixB)

realDataCov<-cov(sample_distribution)
realDataCovB<-cov(sample_distributionB)
covs<-t(realDataCov[upper.tri(realDataCov,diag = T)])
covsB<-t(realDataCovB[upper.tri(realDataCovB,diag = T)])
dataSSs<-data.frame(covs,covsB)
names(dataSSs)<-c("var_X1","cov_X1_X2","var_X2","cov_X1_X3","cov_X2_X3","var_X3","var_X1B","cov_X1_X2B","var_X2B","cov_X1_X3B","cov_X2_X3B","var_X3B")

#Run the abc, currently 1% kept. 
#Neural network regression adjustment.
abcOutputNNet = abc(target=dataSSs, param=paramSims, sumstat=SSsims,tol=0.01,method="neuralnet")

