# TODO: Add comment
# 
# Author: nk3
###############################################################################

library(limma)
library(methods)
library(car)
library(nlme)
library(nortest)
library(vcd)
library(grid)
library(MASS)
library(logistf)
library(lme4)

setwd("C:\\git2\\IMPC_statsWorkingGroup\\PhenStat\\PhenStatSource\\R")
source("classes.R")
#source("classification.R")
source("diagnostictest.R")
source("FETFramework.R")
source("graphsDataset.R")
source("graphsResults.R")
source("LogisticFramework.R")
source("MMFramework.R")
source("PhenList.R")
source("PhenTestResult.R")
source("RRFramework.R")
source("summaryOutput.R")
#source("testDataset.R")
source("TFFramework.R")
source("transformation.R")
source("vectorOutput.R")


setwd("C:\\git2\\IMPC_statsWorkingGroup\\PhenStat\\PhenStatSD")
source("testDataset.R")
source("vectorOutputSD.R")
source("MMFrameworkAssessSD.R")
source("classificationSD.R")
source("LogisticFrameworkSD.R")



#load some data files for logistic regression data

setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
dataset=read.csv("categorical_example 1_aff3_Xray.csv")
levels(dataset$Skull.Shape)
test <- PhenList(dataset,testGenotype="Aff3/Aff3")
test2 <- LRDataset(test, depVariable="Skull.Shape", abnormalValues="Abnormal")

#old method
result2 <- testDataset(phenList=test2,depVariable="Skull.Shape",method="LR")
summaryOutput(result2)

#new method
result2 <- testDataset(phenList=test2,depVariable="Skull.Shape",method="SD_categorical")
summaryOutput(result2)

result2=startLRModel_SD(test2, depVariable="Skull.Shape", outputMessages=TRUE, pThreshold=0.05)
result3 <- finalLRModel_SD(result2, outputMessages=TRUE)
result2 <- testDataset(phenList=test2,depVariable="Skull.Shape",method="SD_categorical")

vectorOutputSD(result3)

debug(vectorOutputSD)
#problem with 

playing1 <- analysisResults(result3)  
playing1$SDmodel.output.summary

playing1 <- analysisResults(result2)  
playing1$SDmodel.output.summary

