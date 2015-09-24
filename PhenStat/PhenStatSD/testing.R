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

setwd("C:\\git2\\IMPC_statsWorkingGroup\\PhenStat\\PhenStatSource\\R")
source("classes.R")
source("classification.R")
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
source("testDataset.R")
source("TFFramework.R")
source("transformation.R")
source("vectorOutput.R")

library(methods)
library(car)
library(nlme)
library(nortest)
library(MASS)
library(logistf)
library(lme4)


setwd("C:\\git2\\IMPC_statsWorkingGroup\\PhenStat\\PhenStatSD")
source("testDataset.R")
source("vectorOutputSD.R")
source("MMFrameworkAssessSD.R")


setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
data_typical=read.csv("MAHN_DEXAdata.csv")
PhenObject=PhenList(data_typical, testGenotype="Mysm1/+", dataset.colname.sex="Gender", dataset.clean=TRUE,dataset.colname.batch="Assay.Date", dataset.colname.genotype="Genotype", dataset.colname.weight="Weight", dataset.values.male="Male", 		dataset.values.female="Female") 
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight")  

Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight", method="SD_continuous")  
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withWeight", method="SD_continuous")  

data2=read.csv("MAHN_DEXAdataNoWeight.csv")

PhenObject=PhenList(data2, testGenotype="Mysm1/+", dataset.colname.sex="Gender", dataset.clean=TRUE,dataset.colname.batch="Assay.Date", dataset.colname.genotype="Genotype", dataset.colname.weight="Weight", dataset.values.male="Male", 		dataset.values.female="Female") 
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight")  
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight", method="SD_continuous")  
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withWeight", method="SD_continuous")  


#To do
#1 check handles one sex dataset appropriately
#2 double check correctly grabs data in the right places from various models
#3 large scale run locally

