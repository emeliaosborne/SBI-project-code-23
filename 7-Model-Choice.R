#Model Choice

#Generate model indices
table_len=40000
equal_1s_2s <- c(rep(1, table_len/2), rep(2, table_len/2))

#Shuffle model indices
modindex <- sample(equal_1s_2s)

sample_size<-10000
sample_meanvector <- c(0,0,0)    
mySigmaSq=1.25
myLambda=1/500

#Function which generates mvn samples and then summary statistics.

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

#Generate reference table.

delta1Sims<-rexp(table_len,1/500)
deltamSims<-rexp(table_len,1/500)
delta2MinusmSims<-rexp(table_len,1/500)

alpha1Sims<-numeric(table_len)
alpha2Sims<-numeric(table_len)

SSsimsX1Var = numeric(table_len)
SSsimsCovX1X2 = numeric(table_len)
SSsimsX2Var= numeric(table_len)
SSsimsCovX1X3 = numeric(table_len)
SSsimsCovX2X3 = numeric(table_len)
SSsimsX3Var= numeric(table_len)
SSsimsX1VarB = numeric(table_len)
SSsimsCovX1X2B = numeric(table_len)
SSsimsX2VarB= numeric(table_len)
SSsimsCovX1X3B = numeric(table_len)
SSsimsCovX2X3B = numeric(table_len)
SSsimsX3VarB= numeric(table_len)
SSsims<-data.frame(SSsimsX1Var,SSsimsCovX1X2,SSsimsX2Var,SSsimsCovX1X3, SSsimsCovX2X3, SSsimsX3Var,SSsimsX1VarB,SSsimsCovX1X2B,SSsimsX2VarB,SSsimsCovX1X3B, SSsimsCovX2X3B, SSsimsX3VarB)


for (i in 1:table_len){
  alpha1Sims[i]=rbeta(1,2,3)
  if(modindex[i]==1){
    alpha2Sims[i]=alpha1Sims[i]
  }
  else{
    alpha2Sims[i]=rbeta(1,2,3)
    while (abs(alpha2Sims[i]-alpha1Sims[i])<0.1) {
      alpha2Sims[i]=rbeta(1,2,3)
    }
  }
  currentDelta1=delta1Sims[i]
  currentDeltam=deltamSims[i]
  currentDelta2Minusm=delta2MinusmSims[i]
  
  currentAlpha1=alpha1Sims[i]
  currentAlpha2=alpha2Sims[i]
  currentCovMatrix=matrix(c(mySigmaSq*(currentDelta1+currentDeltam+currentDelta2Minusm),
                            currentAlpha1*mySigmaSq*(currentDelta1+currentDeltam), 
                            0, 
                            currentAlpha1*mySigmaSq*(currentDelta1+currentDeltam), 
                            (currentAlpha1^2*(currentDelta1+currentDeltam)+(1-currentAlpha1)^2*(currentDelta1+currentDeltam)+currentDelta2Minusm)*mySigmaSq,
                            (1-currentAlpha1)*mySigmaSq*currentDelta1, 
                            0, 
                            (1-currentAlpha1)*mySigmaSq*currentDelta1, 
                            mySigmaSq*(currentDelta1+currentDeltam+currentDelta2Minusm)),
                          nrow=3,ncol=3,byrow = TRUE)
  currentCovMatrixB=matrix(c(mySigmaSq/(2*myLambda)*(exp(2*myLambda*(currentDelta1+currentDeltam+currentDelta2Minusm))-1),
                             currentAlpha2*mySigmaSq/(2*myLambda)*(exp(2*myLambda*(currentDelta1+currentDeltam))-1), 
                             0, 
                             currentAlpha2*mySigmaSq/(2*myLambda)*(exp(2*myLambda*(currentDelta1+currentDeltam))-1), 
                             (currentAlpha2^2*(exp(2*myLambda*(currentDelta1+currentDeltam))-1)+(1-currentAlpha2)^2*(exp(2*myLambda*(currentDelta1+currentDeltam))-1)+exp(2*myLambda*(currentDelta1+currentDeltam+currentDelta2Minusm))-exp(2*myLambda*(currentDelta1+currentDeltam)))*mySigmaSq/(2*myLambda), 
                             (1-currentAlpha2)*mySigmaSq/(2*myLambda)*(exp(2*myLambda*currentDeltam)-1), 
                             0, 
                             (1-currentAlpha2)*mySigmaSq/(2*myLambda)*(exp(2*myLambda*currentDeltam)-1), 
                             mySigmaSq/(2*myLambda)*(exp(2*myLambda*(currentDelta1+currentDeltam+currentDelta2Minusm))-1)),nrow=3,ncol=3,byrow = TRUE)
  
  
  SSsims[i,]<-dataSimulation(sample_size,currentCovMatrix,currentCovMatrixB)
  
}

#create observed data and get summary statistics

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
                                     mySigmaSq*(myDelta1+myDeltam+myDelta2Minusm)),nrow=3,ncol=3,byrow = TRUE)
sample_covariance_matrixB <- matrix(c(mySigmaSq/(2*myLambda)*(exp(2*myLambda*(myDelta1+myDeltam+myDelta2Minusm))-1),
                                      myAlpha*mySigmaSq/(2*myLambda)*(exp(2*myLambda*(myDelta1+myDeltam))-1), 
                                      0, 
                                      myAlpha*mySigmaSq/(2*myLambda)*(exp(2*myLambda*(myDelta1+myDeltam))-1), 
                                      (myAlpha^2*(exp(2*myLambda*(myDelta1+myDeltam))-1)+(1-myAlpha)^2*(exp(2*myLambda*(myDelta1+myDeltam))-1)+exp(2*myLambda*(myDelta1+myDeltam+myDelta2Minusm))-exp(2*myLambda*(myDelta1+myDeltam)))*mySigmaSq/(2*myLambda), 
                                      (1-myAlpha)*mySigmaSq/(2*myLambda)*(exp(2*myLambda*myDeltam)-1), 
                                      0, 
                                      (1-myAlpha)*mySigmaSq/(2*myLambda)*(exp(2*myLambda*myDeltam)-1), 
                                      mySigmaSq/(2*myLambda)*(exp(2*myLambda*(myDelta1+myDeltam+myDelta2Minusm))-1)),nrow=3,ncol=3,byrow = TRUE)

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

#modindex_vec=as.numeric(modindex)

#ABC using rejection at 1% acceptance.
mc_output=postpr(dataSSs,modindex,SSsims,tol=0.01,method='rejection')

#The estimated posteriors.
c(
  summary(mc_output$values)[1]/(summary(mc_output$values)[1]+summary(mc_output$values)[2]),
  summary(mc_output$values)[2]/(summary(mc_output$values)[1]+summary(mc_output$values)[2])
)

#The estimated Bayes factor.
(summary(mc_output$values)[1]/(summary(mc_output$values)[1]+summary(mc_output$values)[2]))/(summary(mc_output$values)[2]/(summary(mc_output$values)[1]+summary(mc_output$values)[2]))

#ABC using neural network classifier, 1% acceptance.
mc_output_NNet=postpr(dataSSs,modindex,SSsims,tol=0.01,method='neuralnet')

#The estimated posteriors.
c(
  mc_output_NNet$pred[1],
  mc_output_NNet$pred[2]
)

#The estimated Bayes factor.
mc_output_NNet$pred[1]/mc_output_NNet$pred[2]
