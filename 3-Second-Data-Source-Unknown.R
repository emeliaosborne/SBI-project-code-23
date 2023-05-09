#Second Source of Data (Unknown).
#This R script contains code to perform ABC for model with 3 populations, 1 admixture event.
#Uses two sources of data, relationship between unknown

library("abc")
library("MASS")

sample_size<-10000
sample_meanvector <- c(0,0,0) 
mySigmaSq=1.25
mySigmaSqB=1.5

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
SSsimsA<-data.frame(SSsimsX1Var,SSsimsCovX1X2,SSsimsX2Var,SSsimsCovX1X3, SSsimsCovX2X3, SSsimsX3Var)
SSsimsB<-data.frame(SSsimsX1VarB,SSsimsCovX1X2B,SSsimsX2VarB,SSsimsCovX1X3B, SSsimsCovX2X3B, SSsimsX3VarB)

delta1Sims<-rexp(simsSize,1/500)
deltamSims<-rexp(simsSize,1/500)
delta2MinusmSims<-rexp(simsSize,1/500)
delta1BSims<-rexp(simsSize,1/500)
deltamBSims<-rexp(simsSize,1/500)
delta2MinusmBSims<-rexp(simsSize,1/500)
alphaSims<-rbeta(simsSize,2,3)

paramSimsA<-data.frame(delta1Sims,deltamSims,delta2MinusmSims,alphaSims)
paramSimsB<-data.frame(delta1BSims,deltamBSims,delta2MinusmBSims,alphaSims)
paramSims<-data.frame(delta1Sims,deltamSims,delta2MinusmSims,delta1BSims,deltamBSims,delta2MinusmBSims,alphaSims)
names(paramSims)<-c("delta1Sims", "deltamSims","delta2MinusmSims","delta1BSims", "deltamBSims","delta2MinusmBSims", "alphaSims")

for (val in 1:simsSize) {
  currentDelta1=delta1Sims[val]
  currentDeltam=deltamSims[val]
  currentDelta2Minusm=delta2MinusmSims[val]
  currentDelta1B=delta1BSims[val]
  currentDeltamB=deltamBSims[val]
  currentDelta2MinusmB=delta2MinusmBSims[val]
  currentAlpha=alphaSims[val]
  currentCovMatrix=matrix(c(mySigmaSq*(currentDelta1+currentDeltam+currentDelta2Minusm),
                            currentAlpha*mySigmaSq*(currentDelta1+currentDeltam), 
                            0, 
                            currentAlpha*mySigmaSq*(currentDelta1+currentDeltam),
                            (currentAlpha^2*(currentDelta1+currentDeltam)+(1-currentAlpha)^2*(currentDelta1+currentDeltam)+currentDelta2Minusm)*mySigmaSq,
                            (1-currentAlpha)*mySigmaSq*currentDelta1, 
                            0, 
                            (1-currentAlpha)*mySigmaSq*currentDelta1, 
                            mySigmaSq*(currentDelta1+currentDeltam+currentDelta2Minusm)), nrow=3,ncol=3,byrow = TRUE)
  
  currentCovMatrixB=matrix(c(mySigmaSqB*(currentDelta1B+currentDeltamB+currentDelta2MinusmB),
                             currentAlpha*mySigmaSqB*(currentDelta1B+currentDeltamB), 
                             0, 
                             currentAlpha*mySigmaSqB*(currentDelta1B+currentDeltamB),
                             (currentAlpha^2*(currentDelta1B+currentDeltamB)+(1-currentAlpha)^2*(currentDelta1B+currentDeltamB)+currentDelta2MinusmB)*mySigmaSqB,
                             (1-currentAlpha)*mySigmaSqB*currentDelta1B, 
                             0, 
                             (1-currentAlpha)*mySigmaSqB*currentDelta1B, 
                             mySigmaSqB*(currentDelta1B+currentDeltamB+currentDelta2MinusmB)), nrow=3,ncol=3,byrow = TRUE)
  
  SSsims[val,]<-dataSimulation(sample_size,currentCovMatrix,currentCovMatrixB)
  SSsimsA[val,]<-dataSimulationN(sample_size,currentCovMatrix)
  SSsimsB[val,]<-dataSimulationN(sample_size,currentCovMatrixB)
}

#Generate 'observed' data.
myDelta1=350
myDeltam=160
myDelta2Minusm=640
myDelta1B=176
myDeltamB=63
myDelta2MinusmB=402
myAlpha=0.25

sample_covariance_matrix <- matrix(c(mySigmaSq*(myDelta1+myDeltam+myDelta2Minusm),
                                     myAlpha*mySigmaSq*(myDelta1+myDeltam), 
                                     0, 
                                     myAlpha*mySigmaSq*(myDelta1+myDeltam), 
                                     (myAlpha^2*(myDelta1+myDeltam)+(1-myAlpha)^2*(myDelta1+myDeltam)+myDelta2Minusm)*mySigmaSq,
                                     (1-myAlpha)*mySigmaSq*myDelta1, 
                                     0, 
                                     (1-myAlpha)*mySigmaSq*myDelta1, 
                                     mySigmaSq*(myDelta1+myDeltam+myDelta2Minusm)),nrow=3,ncol=3,byrow = TRUE)
sample_covariance_matrixB <- matrix(c(mySigmaSqB*(myDelta1B+myDeltamB+myDelta2MinusmB),
                                      myAlpha*mySigmaSqB*(myDelta1B+myDeltamB), 
                                      0, 
                                      myAlpha*mySigmaSqB*(myDelta1B+myDeltamB), 
                                      (myAlpha^2*(myDelta1B+myDeltamB)+(1-myAlpha)^2*(myDelta1B+myDeltamB)+myDelta2MinusmB)*mySigmaSqB,
                                      (1-myAlpha)*mySigmaSqB*myDelta1B, 
                                      0, 
                                      (1-myAlpha)*mySigmaSqB*myDelta1B, 
                                      mySigmaSqB*(myDelta1B+myDeltamB+myDelta2MinusmB)),nrow=3,ncol=3,byrow = TRUE)

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
dataSSsA<-data.frame(covs)
dataSSsB<-data.frame(covsB)
names(dataSSs)<-c("var_X1","cov_X1_X2","var_X2","cov_X1_X3","cov_X2_X3","var_X3","var_X1B","cov_X1_X2B","var_X2B","cov_X1_X3B","cov_X2_X3B","var_X3B")


#Run the abc, currently 1% kept. 
#Neural network regression adjustment.
abcOutputNNet = abc(target=dataSSs, param=paramSims, sumstat=SSsims,tol=0.01,method="neuralnet")

#ABCs for combined alpha
abcOutputANNet = abc(target=dataSSsA, param=paramSimsA, sumstat=SSsimsA,tol=0.01,method="neuralnet")
abcOutputBNNet = abc(target=dataSSsB, param=paramSimsB, sumstat=SSsimsB,tol=0.01,method="neuralnet")

#posterior of alpha given both sources of data.
muA=mean(abcOutputANNet$adj.values[,4])
muB=mean(abcOutputBNNet$adj.values[,4])
sigmaSqA=var(abcOutputANNet$adj.values[,4])
sigmaSqB=var(abcOutputBNNet$adj.values[,4])

muC=(muA+muB)/2
sigmaSqC=(sigmaSqA+sigmaSqB)/4
