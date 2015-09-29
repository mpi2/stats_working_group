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

#debug(testDataset)

setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
data_typical=read.csv("MAHN_DEXAdata.csv")
PhenObject=PhenList(data_typical, testGenotype="Mysm1/+", dataset.colname.sex="Gender", dataset.clean=TRUE,dataset.colname.batch="Assay.Date", dataset.colname.genotype="Genotype", dataset.colname.weight="Weight", dataset.values.male="Male", 		dataset.values.female="Female") 
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight", method="MM")  
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withWeight", method="SD_continuous")  

summary(Modelfeatures@analysisResults$model.SDmodel_output)
vectorOutputSD(Modelfeatures)   # problem not grabbing the right values

undebug(testDataset)


undebug(vectorOutputSD)

#ruled 
#table length is correct
#components are in the same order

Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight", method="SD_continuous")  

#undebug(analysisResults)


data2=read.csv("MAHN_DEXAdataNoWeight.csv")
PhenObject=PhenList(data2, testGenotype="Mysm1/+", dataset.colname.sex="Gender", dataset.clean=TRUE,dataset.colname.batch="Assay.Date", dataset.colname.genotype="Genotype", dataset.colname.weight="Weight", dataset.values.male="Male", 		dataset.values.female="Female") 
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight")  
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight", method="SD_continuous")  
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withWeight", method="SD_continuous")  


setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
data3= read.csv("MAHN_DEXAdataNoWeightNoFemale.csv")
PhenObject=PhenList(data3, testGenotype="Mysm1/+", dataset.colname.sex="Gender", dataset.clean=TRUE,dataset.colname.batch="Assay.Date", dataset.colname.genotype="Genotype", dataset.colname.weight="Weight", dataset.values.male="Male", 		dataset.values.female="Female") 
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight")  
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight", method="SD_continuous")  
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withWeight", method="SD_continuous")  

data4=read.csv("MAHN_DEXAdataNoWeightNoKOFemale.csv")
PhenObject=PhenList(data4, testGenotype="Mysm1/+", dataset.colname.sex="Gender", dataset.clean=TRUE,dataset.colname.batch="Assay.Date", dataset.colname.genotype="Genotype", dataset.colname.weight="Weight", dataset.values.male="Male", 		dataset.values.female="Female") 
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight")  
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight", method="SD_continuous")  
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withWeight", method="SD_continuous")  

#correctly halts the analysis. 

#To do
#1 double check correctly grabs data in the right places from various models
#2 large scale run locally

