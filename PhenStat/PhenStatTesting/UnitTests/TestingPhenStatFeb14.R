# TODO: Testing following restructuring of the dataset checks
# # Author: nk3
###############################################################################
library(limma)
library(methods)
library(car)
library(nlme)
library(nortest)
library(vcd)
library(grid)




#download test datasets
setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
data_typical=read.csv("MAHN_DEXAdata.csv")
data_oneGender=read.csv("MAHN_DEXAdata_1Gender.csv")
data_3Genders=read.csv("MAHN_DEXAdata_3genders.csv")
data_genderas1or2=read.csv("MAHN_DEXAdata_genderas1or2.csv")
data_noWeight=read.csv("MAHN_DEXAdata_noWeight.csv")
data_noAssayDate=read.csv("MAHN_DEXAdata-noAssayDateColumn.csv")
data_noGendercolumn=read.csv("MAHN_DEXAdata-noGenderColumn.csv")
data_3Genotypes=read.csv("MAHN_DEXAdata_3genotypes.csv")


#standard dataframe testing
boxplotGenderGenotype(PhenObject, "Lean.Mass", "Lean Mass (g)")  #appropriate error 
PhenObject=PhenList(data_typical) # appropriate error
PhenObject=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=FALSE)  # good errors
#good still removed empty "" genotype level
PhenObject=PhenList(data_typical, testGenotype="Mysm1/+") #good messages
PhenObject=PhenList(data_3Genotypes, testGenotype="Mysm1/+") #test 3 genotype test dataset

Modelfeatures=testDataset(PhenObject)   #test you haven't defined depVariable
Modelfeatures=testDataset(PhenObject, depVariable="Lean.mass") #test if you have typed the name incorrectly
Modelfeatures=testDataset(PhenObject, depVariable="Lean.Mass")
Modelfeatures=testDataset(PhenObject, depVariable="Lean.Mass", equation="withouWeight") # test entered equation option incorerctly
Modelfeatures=testDataset(PhenObject, depVariable="Lean.mass", equation="withoutWeight") # test entered dep Variable incorerctly
Modelfeatures=testDataset(PhenObject, depVariable="Age.In.Weeks", equation="withWeight") #test that when weight is not significant reverts to equation without weight and reports appropriately.
Modelfeatures
Modelfeatures=testDataset(PhenObject, depVariable="Lean.Mass", equation="withoutWeight")
Modelfeatures
PhenObject=PhenList(data_3Genotypes, testGenotype="Mysm1/+", dataset.clean=FALSE) # impact of turning of clean 
PhenObject=PhenList(data_typical,testGenotype="Mysm1/+", dataset.clean=FALSE) # impact of turning of clean 
Modelfeatures=testDataset(PhenObject, depVariable="Lean.Mass", equation="withoutWeight")  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PROBLEM HERE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Modelfeatures

Modelfeatures=testDataset(PhenObject, depVariable="Lean.Mass", equation="withWeight")
Modelfeatures
######################################
#OneGender Dataset testing
levels(data_oneGender$Genotype)
levels(data_oneGender$Gender)
PhenObject=PhenList(data_oneGender, dataset.clean=TRUE, testGenotype="Mysm1/+") 
Modelfeatures=testDataset(PhenObject, depVariable="Lean.Mass", equation="withWeight") ~~~~~~~~~~~~~~~~~~~~~~~~PROBLEM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PhenObject=PhenList(data_oneGender, testGenotype="Mysm1/+", dataset.clean=FALSE)
Modelfeatures=testDataset(PhenObject, depVariable="Lean.Mass", equation="withWeight")


######################################
#3 gender testing
levels(data_3Genders$Genotype)
PhenObject2=PhenList(data_3Genders, dataset.clean=TRUE, testGenotype="Mysm1/+")
Modelfeatures=testDataset(PhenObject2, depVariable="Lean.Mass", equation="withWeight")
Modelfeatures=testDataset(PhenObject2, depVariable="Lean.Mass", equation="withoutWeight")
PhenObject2=PhenList(data_3Genders, dataset.clean=FALSE, testGenotype="Mysm1/+")

#######################################################
# gender as 1 or 2
levels(data_genderas1or2$Genotype)
levels(data_genderas1or2$Gender)
levels(data_typical$Gender)
data_genderas1or2$Gender
PhenObject2=PhenList(data_genderas1or2, dataset.clean=TRUE, testGenotype="Mysm1/+")  # nice message
#works nicely as in error clear need to relabel male and female
PhenObject2=PhenList(data_genderas1or2, dataset.clean=FALSE, testGenotype="Mysm1/+")
PhenObject=PhenList(data_genderas1or2, dataset.clean=TRUE, dataset.values.female="1", dataset.values.male="2", testGenotype="Mysm1/+")
Modelfeatures=testDataset(PhenObject2, depVariable="Lean.Mass", equation="withWeight")
Modelfeatures=testDataset(PhenObject2, depVariable="Lean.Mass", equation="withoutWeight")


PhenObject2=PhenList(data_genderas1or2, dataset.clean=FALSE, dataset.values.female="1", testGenotype="Mysm1/+")
PhenObject2=PhenList(data_genderas1or2, dataset.clean=FALSE, dataset.values.female=1, dataset.values.male=2, testGenotype="Mysm1/+")
PhenObject2=PhenList(data_genderas1or2, dataset.clean=TRUE, dataset.values.female=1, dataset.values.male=2, testGenotype="Mysm1/+")
#need to clearly give example of this in code as you need datset.clean=TRUE as you are modifying dataset. ##DONE

data_genderFunnyStrings=read.csv("MAHN_DEXAdata_genderasFunnyNamesAsStrings.csv")
PhenObject2=PhenList(data_genderFunnyStrings, dataset.clean=TRUE, testGenotype="Mysm1/+")
levels(data_genderFunnyStrings$Gender)
PhenObject=PhenList(data_genderFunnyStrings, dataset.clean=TRUE, dataset.values.female="FemaleHonest", dataset.values.male="MaleHonest", testGenotype="Mysm1/+")
Modelfeatures=testDataset(PhenObject, depVariable="Lean.Mass", equation="withWeight")
Modelfeatures

################################################
#no weight
data_noWeight=read.csv("MAHN_DEXAdata_noWeight.csv")
levels(data_noWeight$Genotype)
PhenObject2=PhenList(data_noWeight, dataset.clean=TRUE, testGenotype="Mysm1/+")
Modelfeatures=testDataset(PhenObject2, depVariable="Lean.Mass", equation="withWeight") # good error that had to use withoutWeight  ~~~~~~~~~~~~~~~~~~~problem odd warning message now appeared~~~~~~~~~~~~~~~~~~~~~~~
Modelfeatures

PhenObject2=PhenList(data_noWeight, dataset.clean=FALSE, testGenotype="Mysm1/+")
Modelfeatures=testDataset(PhenObject2, depVariable="Lean.Mass", equation="withoutWeight")
Modelfeatures

PhenObject2=PhenList(data_noWeight, dataset.clean=TRUE,dataset.colname.weight="Est.Total.Tissue.Mass", testGenotype="Mysm1/+")  ~~~~~~~~~~~~~~~~~problem odd error message~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Modelfeatures=testDataset(PhenObject2, depVariable="Lean.Mass", equation="withWeight") 
Modelfeatures

##########################################
#no Assay date
PhenObject=PhenList(data_noAssayDate, dataset.clean=TRUE, testGenotype="Mysm1/+")
Modelfeatures=testDataset(PhenObject, depVariable="Lean.Mass", equation="withoutWeight")
Modelfeatures

#################################################
#no gender column
PhenObject=PhenList(data_noGendercolumn, dataset.clean=TRUE, testGenotype="Mysm1/+")		

#############################
PhenList(dataset, testGenotype, refGenotype="+/+", hemiGenotype=NULL, 
		outputMessages=TRUE, dataset.clean=TRUE, 
		dataset.colname.batch=NULL, dataset.colname.genotype=NULL, 
		dataset.colname.gender=NULL, dataset.colname.weight=NULL, 
		dataset.values.missingValue=NULL, dataset.values.male=NULL, 
		dataset.values.female=NULL, dataset.stat=NULL)


PhenObject=PhenList(data_typical,testGenotype="Mysm1/+",dataset.stat=TRUE ) 
PhenObject





#########################################
#test graphsDataset
#########################################
#1: boxplotGenderGenotype
############################
#test standard dataset
PhenObject=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=FALSE)
boxplotGenderGenotype(PhenObject, "Lean.Mass", "Lean Mass (g)")
boxplotGenderGenotype(PhenObject, "Lean.mass", "Lean Mass (g)")  # good error message

#test one gender dataset
PhenObject=PhenList(data_oneGender, dataset.clean=TRUE, testGenotype="Mysm1/+") 
boxplotGenderGenotype(PhenObject, "Lean.Mass", "Lean Mass (g)")
boxplotGenderGenotype(PhenObject,  "Lean Mass (g)")
boxplotGenderGenotype(PhenObject, graphingName="Lean Mass (g)") #test omitting the PhenObject
boxplotGenderGenotype(PhenObject, "Lean.Mass")

#2: boxplotGenderGenotypeBatch
###############################
PhenObject=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=FALSE)
PhenObject=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=TRUE)
boxplotGenderGenotypeBatch(PhenObject, "Lean.Mass", "Lean Mass (g)")
boxplotGenderGenotypeBatch(PhenObject, "Lean.mass", "Lean Mass (g)")

#test one gender dataset
PhenObject=PhenList(data_oneGender, dataset.clean=TRUE, testGenotype="Mysm1/+") 
boxplotGenderGenotypeBatch(PhenObject, "Lean.Mass", "Lean Mass (g)")
boxplotGenderGenotypeBatch(PhenObject,  "Lean Mass (g)")
boxplotGenderGenotypeBatch(PhenObject, graphingName="Lean Mass (g)") #test omitting the PhenObject
boxplotGenderGenotypeBatch(PhenObject, "Lean.Mass")

#3 scatterplotGenotypeWeight
############################
PhenObject=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=FALSE)
scatterplotGenotypeWeight(PhenObject, "Lean.Mass", "Lean Mass (g)")

PhenObject2=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=TRUE)
scatterplotGenotypeWeight(PhenObject2, "Lean.Mass", "Lean Mass (g)")
scatterplotGenotypeWeight(PhenObject2, "Lean.mass", "Lean Mass (g)")
#good error that depVariable column is not in dataset
scatterplotGenotypeWeight(PhenObject,  "Lean Mass (g)")
scatterplotGenotypeWeight(PhenObject, graphingName="Lean Mass (g)") #test omitting the PhenObject
scatterplotGenotypeWeight(PhenObject, "Lean.Mass")

PhenObject=PhenList(data_oneGender, dataset.clean=TRUE, testGenotype="Mysm1/+") 
scatterplotGenotypeWeight(PhenObject, "Lean.Mass", "Lean Mass (g)")


#what about when lots of missing value in weight variable?

setwd("X:\\2013\\BioconductorPackageTesting\\BugTesting\\11thFeb14")
data1=read.csv("Error in ctrfn(levels(x), contrasts = contrasts)contrasts not defined for 0 degrees of freedom.csv")
colnames(data1)
PhenObject2=PhenList(data1, testGenotype="3", refGenotype="1", dataset.clean=TRUE, dataset.colname.weight="Weight")
scatterplotGenotypeWeight(PhenObject2, "Param", "Param")  #good error message
Modelfeatures=testDataset(PhenObject2, depVariable="Param", equation="withWeight")
Modelfeatures=testDataset(PhenObject2, depVariable="Param", equation="withoutWeight")

data2=read.csv("Error in cvm.test(Gp2$res)  sample size must be greater than 7.csv")
PhenObject2=PhenList(data2, testGenotype="3", refGenotype="1", dataset.clean=TRUE, dataset.colname.weight="Weight")
scatterplotGenotypeWeight(PhenObject2, "Param", "Param")  #good error message
Modelfeatures=testDataset(PhenObject2, depVariable="Param", equation="withWeight")
Modelfeatures=testDataset(PhenObject2, depVariable="Param", equation="withoutWeight")

####################################################
#graphsResults.R
########################################################
#1: qqplotGenotype<-function(phenList, phenTestResult)
PhenObject2=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=TRUE)
result <- testDataset(PhenObject2,depVariable="Bone.Area")
qqplotGenotype(phenTestResult=result)
result <- testDataset(PhenObject2,depVariable="Lean.Mass")
qqplotGenotype( phenTestResult=result)
#testing a one gender dataset
PhenObject=PhenList(data_oneGender, dataset.clean=TRUE, testGenotype="Mysm1/+") 
result <- testDataset(PhenObject,depVariable="Lean.Mass")
qqplotGenotype(phenTestResult=result) # not tested as a one gender dataset and problem higher up with one gneder dataset

#2 qqplotRandomEffects<-function(phenList, phenTestResult, keep_batch=NULL)
PhenObject2=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=TRUE)
result <- testDataset(PhenObject2,depVariable="Bone.Area")
qqplotRandomEffects(phenTestResult=result)

#testing a dataset where batch is not releveant
result <- testDataset(PhenObject2,depVariable="Est.Total.Tissue.Mass")
result$model.effect.batch
qqplotRandomEffects(phenTestResult=result) # good diagnostics

#testing a one gender dataset
PhenObject=PhenList(data_oneGender, testGenotype="Mysm1/+", dataset.clean=TRUE)
result <- testDataset(PhenObject,depVariable="Lean.Mass")
qqplotRandomEffects(phenList=PhenObject2, phenTestResult=result) # can't test yet as have problems with one gender dataset
qqplotRandomEffects(phenTestResult=result)
# fine with oneGender datasets

#3 qqplotRotatedResiduals<-function(phenList, phenTestResult, keep_batch=NULL)
PhenObject2=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=TRUE)
result <- testDataset(PhenObject2,depVariable="Bone.Area")
qqplotRotatedResiduals(phenTestResult=result)
result <- testDataset(PhenObject2,depVariable="Est.Total.Tissue.Mass")
qqplotRotatedResiduals(phenTestResult=result)  # not plotted as batch not significant
#testing a one gender dataset
PhenObject=PhenList(data_oneGender, testGenotype="Mysm1/+", dataset.clean=TRUE)
result <- testDataset(PhenObject,depVariable="Lean.Mass")
qqplotRotatedResiduals(phenTestResult=result) #can't test yet as problem with one gender datasets


#4 boxplotResidualBatch<-function(phenList, phenTestResult)
PhenObject2=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=TRUE)
result <- testDataset(PhenObject2,depVariable="Bone.Area")
result <- testDataset(PhenObject2,depVariable="Bone.Area", equation="withoutWeight")
boxplotResidualBatch(phenTestResult=result)

result <- testDataset(PhenObject2,depVariable="Est.Total.Tissue.Mass")
boxplotResidualBatch(phenTestResult=result)

#testing a one gender dataset
PhenObject=PhenList(data_oneGender, testGenotype="Mysm1/+", dataset.clean=TRUE)
result <- testDataset(PhenObject,depVariable="Lean.Mass")
boxplotResidualBatch(phenList=PhenObject, phenTestResult=result) # works fine with a one gender dataset  #can't test due to rpblem with one gender datasets


#5 plotResidualPredicted<-function(phenList, phenTestResult)
PhenObject2=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=TRUE)
result <- testDataset(PhenObject2,depVariable="Bone.Area")
result <- testDataset(PhenObject2,depVariable="Bone.Area", equation="withoutWeight")
plotResidualPredicted(phenTestResult=result)
result <- testDataset(PhenObject2,depVariable="Est.Total.Tissue.Mass")
plotResidualPredicted(phenList=PhenObject2, phenTestResult=result)
#testing a one gender dataset
PhenObject=PhenList(data_oneGender, testGenotype="Mysm1/+", dataset.clean=TRUE)
result <- testDataset(PhenObject,depVariable="Lean.Mass")
plotResidualPredicted(phenTestResult=result) # can't test

####################################################
#classification.R= classificationTag
####################################################
#classificationTag<-function(result, interactionMode=TRUE, phenotypeThreshold=0.01)
setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
data_typical=read.csv("MAHN_DEXAdata.csv")
data_oneGender=read.csv("MAHN_DEXAdata_1Gender.csv")

PhenObject2=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=TRUE)
result <- testDataset(PhenObject2,depVariable="Lean.Mass")
vectorOutput(result)
result
#genotype effect: 0.0006 - both genders equally
classificationTag(PhenObject2) #check if input the wrong sort of item  -fails
classificationTag(result, userMode=TRUE, phenotypeThreshold=0.01)#checks if input the right options for userMode
classificationTag(result, userMode="summaryOutput", phenotypeThreshold=0.001) #"With phenotype threshold value 0.001 - no significant change"
classificationTag(result, userMode="vectorOutput")# "If phenotype is significant it is for the one genotype tested"
classificationTag(result, userMode="SummaryOutput", phenotypeThreshold=0.05) #check error if type wrongly
classificationTag(result, userMode="summaryOutput", phenotypeThreshold=0.7)#With phenotype threshold value 0.7 - no significant change"

PhenObject2=PhenList(data_oneGender, testGenotype="Mysm1/+", dataset.clean=TRUE)
result2 <- testDataset(PhenObject2,depVariable="Lean.Mass")
classificationTag(result2, userMode="summaryOutput", phenotypeThreshold=0.01) #correct "With phenotype threshold value 0.01 - no significant change" #########CANT CHECK

classificationTag(result2, userMode="vectorOutput")  #correct "If phenotype is significant it is for the one genotype tested" ########CANT CHECK

PhenObject2=PhenList(data_oneGender, testGenotype="Mysm1/+", dataset.clean=TRUE)
result2 <- testDataset(PhenObject2,depVariable="Bone.Mineral.Content")
classificationTag(result2, userMode="summaryOutput", phenotypeThreshold=0.05) #correct "With phenotype threshold value 0.05 - a significant change for the one genotype tested"
classificationTag(result2, userMode="summaryOutput", phenotypeThreshold=0.01) #correct With phenotype threshold value 0.01 - no significant change"
classificationTag(result2, userMode="vectorOutput")  #correct "If phenotype is significant it is for the one genotype tested"


####################################################
#vectorOutput.R vectorOutput(phenTestResult)
####################################################
setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
data_typical=read.csv("MAHN_DEXAdata.csv")
data_oneGender=read.csv("MAHN_DEXAdata_1Gender.csv")
PhenObject2=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=TRUE)
result <- testDataset(PhenObject2,depVariable="Lean.Mass")
vectorOutput(result) 

PhenObject=PhenList(data_oneGender, testGenotype="Mysm1/+", dataset.clean=TRUE) ########CANT CHECK
result <- testDataset(PhenObject,depVariable="Lean.Mass")
vectorOutput(result) 
result


####################################################
#vsummaryOutput
####################################################
#1. summaryOutput <- function(phenTestResult,phenotypeThreshold=0.01)
setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
data_typical=read.csv("MAHN_DEXAdata.csv")
data_oneGender=read.csv("MAHN_DEXAdata_1Gender.csv")
PhenObject2=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=TRUE)
result <- testDataset(PhenObject2,depVariable="Lean.Mass")
summaryOutput(result) 

PhenObject=PhenList(data_oneGender, testGenotype="Mysm1/+", dataset.clean=TRUE)  ########CANT CHECK
result <- testDataset(PhenObject,depVariable="Lean.Mass")
summaryOutput(result) 


########################
setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
data_typical=read.csv("MAHN_DEXAdata.csv")
PhenObject=PhenList(data_typical, testGenotype="Mysm1/+", dataset.clean=TRUE)
output2=startModel(phenList=PhenObject, depVariable="Lean.Mass", equation="withoutWeight", outputMessages=FALSE) 
model1=finalModel(phenTestResult=output2, outputMessages = TRUE)
summaryOutput(model1)

boxplotGenderGenotypeBatch(phenList=PhenObject,depVariable="Lean.Mass",graphingName="Lean Mass")
output2=startModel(phenList=PhenObject, equation = "withoutWeight", depVariable="Lean.Mass", keepList=c(keepBatch=TRUE,keepVariance=TRUE,keepWeight=FALSE,keepGender=TRUE,keepInteraction=TRUE))
#note over-wrote the keep Interaction so which ones can you truely fix?
model1=finalModel(phenTestResult=output2, outputMessages = TRUE)
summaryOutput(model1)

##########################
#testing categorical tools
#############################
setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
test <- PhenList(dataset=read.csv("categorical_example 1_aff3_Xray.csv"),testGenotype="Aff3/Aff3")
result <- testDataset(test,depVariable="Skull.Shape",method="FE")
categoricalBarplot(result)
classificationTag(result)
summaryOutput(result)
vectorOutputMatrices(result)
vectorOutput(result)
dataset=read.csv("categorical_example 1_aff3_Xray.csv")
result <- testDataset(test,depVariable="Number.Of.Ribs.Right",method="MM")  # good a stop insufficient variation
result <- testDataset(test,depVariable="Number.Of.Ribs.Right",method="FE") 
categoricalBarplot(result)
summaryOutput(result)
vectorOutput(result)
vectorOutputMatrices(result)

dataset$Tibia.Length

result <- testDataset(test,depVariable="Tibia.Length",method="FE") # good error now that there is not enough data
result <- testDataset(test,depVariable="Tibia.Length",method="MM")  #~~~~~~~~~~~~~~~~~~~~~~~~problem now~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

result <- testDataset(test,depVariable="Number.Of.Caudal.Vertebrae",method="FE")
summaryOutput(result)
vectorOutput(result)
categoricalBarplot(result)
vectorOutputMatrices(result)

test <- PhenList(dataset=read.csv("categorical_example 1_aff3_Xray.csv"),testGenotype="Aff3/Aff3")
result <- testDataset(test,depVariable="Number.Of.Caudal.Vertebrae",method="MM")  
boxplotGenderGenotypeBatch(test, "Number.Of.Caudal.Vertebrae", "Number.Of.Caudal.Vertebrae")
boxplotGenderGenotypeBatch(result, "Number.Of.Caudal.Vertebrae", "Number.Of.Caudal.Vertebrae")


boxplotGenderGenotype(result, "Number.Of.Caudal.Vertebrae", "Number.Of.Caudal.Vertebrae")
boxplotResidualBatch(test)
qqplotRotatedResiduals(test)

##one gender dataset########
setwd("Y:\\2013\\BioconductorPackageTesting\\Dataset")
test <- PhenList(dataset=read.csv("categorical_example 1_aff3_Xray_1gender.csv"),testGenotype="Aff3/Aff3")
result <- testDataset(test,depVariable="Skull.Shape",method="FE")
result
categoricalBarplot(result)
classificationTag(result)
vectorOutputMatrices(result)
summaryOutput(result)

result <- testDataset(test,depVariable="Number.Of.Caudal.Vertebrae",method="FE")
result
classificationTag(result)   # problem message incorrect as only one gender

result <- testDataset(test,depVariable="Tibia.Length",method="FE") # good error now that there is not enough data
result <- testDataset(test,depVariable="Tibia.Length",method="MM")  #good error not numeric and then not enough data

result <- testDataset(test,depVariable="Number.Of.Caudal.Vertebrae",method="FE")

summaryOutput(result)
vectorOutput(result)
categoricalBarplot(result)
vectorOutputMatrices(result)
result

result <- testDataset(test,depVariable="Number.Of.Caudal.Vertebrae",method="MM")
vectorOutputMatrices(result)



