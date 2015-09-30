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
source("classificationSD.R")
#debug(testDataset)

setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
data_typical=read.csv("MAHN_DEXAdata.csv")
PhenObject=PhenList(data_typical, testGenotype="Mysm1/+", dataset.colname.sex="Gender", dataset.clean=TRUE,dataset.colname.batch="Assay.Date", dataset.colname.genotype="Genotype", dataset.colname.weight="Weight", dataset.values.male="Male", 		dataset.values.female="Female") 
#Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight", method="MM")  
Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withWeight", method="SD_continuous")  

summary(Modelfeatures@analysisResults$model.SDmodel_output)
vectorOutputSD(Modelfeatures)   # 


Modelfeatures=testDataset(PhenObject, depVariable="Bone.Mineral.Density", equation="withoutWeight", method="SD_continuous")  

#undebug(analysisResults)
debug(classificationTagSD)
classificationTag(phenTestResult=Modelfeatures,	userMode="vectorOutput", outputMessages=TRUE)
classificationTag(phenTestResult=Modelfeatures,	userMode="summaryOutput", outputMessages=TRUE)




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


####large scale testing of the code

setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
data_typical=read.csv("MAHN_DEXAdata.csv")
colnames(data_typical)
levels(data_typical$Genotype)

#set the variables of interest
VariablesOfInterest=c( "Weight" , "Nose.To.Tail.Base.Length" ,"Bone.Mineral.Density"  ,"Bone.Mineral.Content" ,  "Bone.Area" ,"Lean.Mass","Fat.Mass", "Est.Total.Tissue.Mass","Fat.Percentage.Estimate" )
		  
countDataPoints<-function(dataset, depVariable){
	print(sum(is.finite(dataset[ , depVariable])))
}

SDPipeline2Analysis<-function(df, GenotypeInterest, GenderColumnName , GenotypeColumnName, refGenotypeString, BatchName, Colony.Prefix, GeneString, VariableInterest, analysisMethod="SD_continuous" ){
	test=PhenList(df, testGenotype=GenotypeInterest, refGenotype=refGenotypeString,dataset.colname.sex=GenderColumnName,dataset.values.male="Male", 	dataset.values.female="Female",  dataset.clean=TRUE, dataset.colname.batch=BatchName,  dataset.colname.genotype=GenotypeColumnName)
	result <- testDataset(test,depVariable=VariableInterest, method=analysisMethod, dataPointsThreshold=4)
	output=vectorOutputSD(result)
	output=c(output, GenotypeInterest, Colony.Prefix, GeneString)
	return(output)
}


#SDPipeline2Analysis(df=data_typical, GenotypeInterest="Mysm1/+",GenderColumnName="Gender", refGenotypeString="+/+", GenotypeColumnName="Genotype", BatchName="Assay.Date", Colony.Prefix="MAHN", GeneString="Mysm1/+", VariableInterest="Lean.Mass", analysisMethod="SD_continuous" )


PhenStatAcrossVariables<-function(data, OutputName,analysisMethod, variablesToTest, ColonyPrefix, BatchName,  GenderColumnName ,GenotypeColumnName,refGenotypeString, GenotypeInterest, GeneName){
	
	#empty object to colleate results too
	finalOutput=c()
		
	for (bob in variablesToTest){
		print(bob)
		
					tryCatch( 
					{PS_output=SDPipeline2Analysis(df=data, GenotypeInterest=as.character(GenotypeInterest),GenderColumnName, GenotypeColumnName ,refGenotypeString, BatchName, Colony.Prefix=ColonyPrefix, GeneString=GeneName, VariableInterest=bob, analysisMethod)
						print(PS_output)
						PS_output=c(PS_output, as.character(GeneName), as.character(ColonyPrefix), as.character(GenotypeInterest))
						
						finalOutput=rbind(PS_output, finalOutput)
					},					
					
					error=function(e){
						
						filename=paste(paste(bob, GenotypeInterest, analysisMethod, bob,"Failuremethod",ColonyPrefix,  sep="_"), "csv", sep=".")
						write.csv(x=PS_output, file=filename)
					})	
			
		}
		write.csv(finalOutput, paste(OutputName, "csv", sep="."))	
	}
				


PhenStatAcrossVariables(data=data_typical, OutputName="testingSDpipeline2",analysisMethod="SD_continuous", variablesToTest=VariablesOfInterest, ColonyPrefix="MAHN", BatchName="Assay.Date",  GenderColumnName="Gender" ,GenotypeColumnName="Genotype",refGenotypeString="+/+", GenotypeInterest="Mysm1/+", GeneName="Mysm1/+")

data2=read.csv("Akt2dataset.csv")
colnames(data2)
levels(data2$Genotype)
PhenStatAcrossVariables(data=data2, OutputName="testingSDpipeline2",analysisMethod="SD_continuous", variablesToTest=VariablesOfInterest, ColonyPrefix="MAAE", BatchName="Assay.Date",  GenderColumnName="Gender" ,GenotypeColumnName="Genotype",refGenotypeString="+/+", GenotypeInterest="Akt2/Akt2", GeneName="Akt2/Akt2")
