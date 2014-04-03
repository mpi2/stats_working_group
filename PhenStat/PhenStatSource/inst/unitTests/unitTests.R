# testing PhenStat 
# Author: nk3, natacourby
##########################################################################
####################General Testing#######################################

# Test data
# UNCOMMENT WHEN NEEDED!!!!!#
# setwd("X:\\2013\\BioconductorPackageTesting\\Dataset")
###############################################################################
# basic checks
test_columnChecks <- function() {    
    data_categorical=read.csv("PhenStat/inst/extdata/test_categorical.csv")
    test <- PhenList(dataset=data_categorical,testGenotype="Aff3/Aff3")
    output <- columnChecks(test$dataset, "Weight", 4) 
    output2 <- columnChecks(test$dataset, "Number.Of.Digits", 8) 
    checkEquals(length(output),3, "check column function returns incorrect output")
    checkEquals(length(output2),3, "check column function returns incorrect output")
    checkTrue(!output[1], "check column function returns incorrect output - there are no column Weight")
    checkTrue(!output[2], "check column function returns incorrect output - there are no column Weight")
    checkTrue(!output[3], "check column function returns incorrect output - there are no column Weight")
    
    checkTrue(output2[1], "check column function returns incorrect output - there is column Number.Of.Digits")
    checkTrue(output2[2], paste("check column function returns incorrect output ",
                    "- column Number.Of.Digits contains numerical values",sep=""))
    checkTrue(!output[3], paste("check column function returns incorrect output ",
                    "- column Number.Of.Digits does not have enough datapoints for genotype/sex combinations",sep=""))
}

###############################################################################
# PhenList object
test_PhenList <- function() {    
    checkException(PhenList(),"PhenList fails with no inputs without exception")
}

###############################################################################
# PhenList with one sex
test_oneSex <- function(){
    data_oneSex=read.csv("PhenStat/inst/extdata/test1_1sex.csv")
    test <- PhenList(dataset=data_oneSex,testGenotype="Mysm1/+")
    checkEquals(length(levels(test$dataset$Sex)),1,paste("length of sex levels is not equals to 1",
                    " - it realises there is only one gender",sep=""))
}

###############################################################################
# Fisher Exact Test generates appropriate classificationTag output
test_FEOutput <- function(){
    data_categorical=read.csv("PhenStat/inst/extdata/test_categorical.csv")
    test <- PhenList(dataset=data_categorical,testGenotype="Aff3/Aff3")
    result <- testDataset(test,depVariable="Skull.Shape",method="FE")
    checkEquals(length(classificationTag(result)),1, "length of classificationTag output is not 1")
    checkTrue(is.character(classificationTag(result)), "classificationTag output is not character")
    checkEquals(classificationTag(result),
            "With phenotype threshold value 0.01 - significant in males and in combined dataset",
            "classificationTag output is incorrect for this dataset")
}

###############################################################################
# testDataset functionality - exception when there are not enough data, 
# switch to Fisher Exact Test when there are categorical data to analyse
test_testDataset <- function(){
    data_categorical=read.csv("PhenStat/inst/extdata/test_categorical.csv")
    test <- PhenList(dataset=data_categorical,testGenotype="Aff3/Aff3")
    
    # Switch from MM to FE for categorical data
    checkEquals(testDataset(test,depVariable="Skull.Shape",method="MM")$method, "FE", 
            "testDataset fails to switch to FE from MM for categorical data analysis")
    
    # All datapoints NA - check for MM
    checkException(testDataset(test,depVariable="Tibia.Length",method="MM"), 
            "testDataset fails to check for sufficient data in MM framework")
    
    # All datapoints NA - check for FE
    checkException(testDataset(test,depVariable="Tibia.Length",method="FE"),
            "testDataset fails with check for sufficient data in FE framework")
}

###############################################################################
# vectorOutput
test_vectorOutput <- function(){
    #file <- system.file("extdata", "test4.csv", package="PhenStat")
    data_typical=read.csv("PhenStat/inst/extdata/test4.csv")
    test <- PhenList(dataset=data_typical,testGenotype="Mysm1/+", dataset.clean=TRUE)
    result <- testDataset(test,depVariable="Lean.Mass")
    vector_results <- vectorOutput(result)   
    checkEquals(length(vector_results), 32, "length of vector is not 32") 
    checkTrue(!is.null(names(vector_results)),"names of vector are missed")
    checkTrue(all(lapply(vector_results,is.character)==TRUE),"not all vector elements are of type character")
}

