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

setwd("U:\\git\\stats_working_group3\\PhenStat\\PhenStatSource\\R")
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


setwd("U:\\git\\stats_working_group3\\PhenStat\\PhenStatSD")
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


###large scale testing
###Load data into R for analysis
setwd("Y:\\Natasha\\Categorical variables\\Logistic Regression - assessing genotype effects\\large scale testing")
data1=read.csv("Oct29_2014_Xray_B6N_MGPSelect_Oor1.csv")

#build an index of KO to drive analysis
Index1=read.csv("Index.csv")

library(plyr)

countDataPoints<-function(dataset, depVariable){
	print(sum(is.finite(dataset[ , depVariable])))
}


XrayVariableList=c("Number.Of.Thoracic.Vertebrae","Number.Of.Lumbar.Vertebrae"  , "Number.Of.Pelvic.Vertebrae","Number.Of.Caudal.Vertebrae"    
		, "Transitional.Vertebrae" , "Shape.Of.Vertebrae"  , "Fusion.Of.Vertebrae", "Processes.On.Vertebrae"  , "Maxilla"  , "Zygomatic.Bone", "Number.Of.Cervical.Vertebrae","Skull.Shape"        
		, "Number.Of.Ribs.Right", "Number.Of.Ribs.Left"  , "Shape.Of.Ribcage" , "Shape.Of.Ribs"    , "Rib.Fusions"  , "Clavicle"  , "Scapula"  ,"Humerus"                       
		, "Radius"  , "Ulna"    ,"Pelvis"   , "Femur"  , "Tibia" , "Fibula"   , "Joints" , "Shape.Of.Spine"    , "Teeth", "Mandible"                     
		,"Number.Of.Digits" ,"Digit.Integrity"  , "Syndactylism" , "Polysyndactylism"   ,"Brachydactylism" , "Kyphosis" , "Lordosis" ,"Scoliosis", "Spinous.Processes" , "Transverse.Processes"  , "Fusion.Processes"              
		, "Caudal.Processes"   , "Cervical.Processes"   , "Lumbar.Processes"   , "Sacral.Processes"  ,"Thoracic.Processes"  ) 		


data1=read.csv("Oct29_2014_Xray_B6N_MGPSelect_Oor1.csv")

#PhenStat wrapper to generate output


LR_Analysis<-function(df, GenotypeInterest, GenotypeName , BatchName, Colony.Prefix,GeneString, VariableInterest ){
	test=PhenList(df, testGenotype=GenotypeInterest, refGenotype="WT",dataset.colname.sex="Gender",dataset.values.male="Male", 
			dataset.values.female="Female",  dataset.clean=TRUE, dataset.colname.batch=BatchName,  dataset.colname.genotype=GenotypeName)
	result <- testDataset(test,depVariable=VariableInterest, method="SD_categorical")
	output=vectorOutputSD(result)
	output=c(output, GeneString, Colony.Prefix)
	
	return(output)
	
}

WT=subset(data1, data1$Genotype2=="WT")
KO=subset(data1, data1$Genotype2!="WT")

LR_Acrosslines<-function(controlData, knockoutData, IndexFile, OutputName){
	
	#empty object to colleate results too
	LR_Output=c()
	
	IndexLength=length(IndexFile[ ,2])
	TestsToRun=c(1:IndexLength)
	
	for (bob in TestsToRun){
		
		#identify the control mice and select them
		ColonyPrefix=IndexFile[bob, 2]
		GeneName=IndexFile[bob, 3]
		Bob_zygosity=IndexFile[bob, 5]
		
		print(ColonyPrefix)
		print(GeneName)
		print(Bob_zygosity)		
		
		#assemble knockout data based on Index criteria
		
		KO_Data_bob=subset(knockoutData, knockoutData[ ,"Colony.Prefix"]==as.character(ColonyPrefix) & knockoutData[ ,"GeneName2"]==as.character(GeneName)& knockoutData[ ,"Genotype2"]==as.character(Bob_zygosity))
		
		#bind - but need same columns
		dataForLR=rbind(KO_Data_bob, controlData)
		
		for(apple in XrayVariableList){
			#using zygosity as genotype for study
			#run the LR analysis
			tryCatch(
					{output=LR_Analysis(dataForLR, GenotypeInterest=as.character(Bob_zygosity), GenotypeName="Genotype2" , BatchName="Assay.Date", Colony.Prefix=ColonyPrefix, GeneString=GeneName,VariableInterest=apple )
						
						print(ColonyPrefix)
						print(Bob_zygosity)
						print(output)
						
						LR_Output=rbind(LR_Output, output)},					
					
					error=function(e){
						
						filename=paste(paste(bob,"FailureLRmethod",Bob_zygosity,ColonyPrefix,  sep="_"), "csv", sep=".")
						write.csv(x=LR_Output, file=filename)
					})	
			
		}
		
	}
	write.csv(LR_Output, paste(OutputName, "csv", sep="."))				
}


debug(LR_Acrosslines)
undebug(LR_Analysis)
LR_Acrosslines(controlData=WT, knockoutData=KO, IndexFile=Index1, OutputName="Playing1.csv")


