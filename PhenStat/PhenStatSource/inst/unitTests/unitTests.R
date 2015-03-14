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
    data_categorical <- read.csv("PhenStat/inst/extdata/test_categorical.csv")
    test <- PhenList(dataset=data_categorical,testGenotype="Aff3/Aff3",outputMessages=FALSE)
    output <- columnChecks(dataset(test), "Weight", 4) 
    output2 <- columnChecks(dataset(test), "Number.Of.Digits", 8) 
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
# PhenList object with one sex
test_oneSex <- function(){
    data_oneSex <- read.csv("PhenStat/inst/extdata/test1_1sex.csv")
    test <- PhenList(dataset=data_oneSex,testGenotype="Mysm1/+",outputMessages=FALSE)
    checkEquals(length(levels(dataset(test)$Sex)),1,paste("length of sex levels is not equals to 1",
                    " - it realises there is only one gender",sep=""))
}
###############################################################################
# Fisher Exact Test generates appropriate classificationTag output
test_FEOutput <- function(){
    data_categorical <- read.csv("PhenStat/inst/extdata/test_categorical.csv")
    test <- PhenList(dataset=data_categorical,testGenotype="Aff3/Aff3",outputMessages=FALSE)
    result <- testDataset(test,depVariable="Skull.Shape",method="FE",outputMessages=FALSE)
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
    test <- PhenList(dataset=data_categorical,testGenotype="Aff3/Aff3",outputMessages=FALSE)    
    # Switch from MM to FE for categorical data
    checkEquals(method(testDataset(test,depVariable="Skull.Shape",method="MM")), "FE", 
            "testDataset fails to switch to FE from MM for categorical data analysis")    
    # All datapoints NA - check for MM
    checkException(testDataset(test,depVariable="Tibia.Length",method="MM"), 
            "testDataset fails to check for sufficient data in MM framework")    
    # All datapoints NA - check for FE
    checkException(testDataset(test,depVariable="Tibia.Length",method="FE"),
            "testDataset fails with check for sufficient data in FE framework")    
    #Not enough data in control subset (requires at least 60)
    checkException(testDataset(test,depVariable="Tibia.Length",method="RR"),
            "testDataset fails with check for sufficient data in RR framework")
}
###############################################################################
# vectorOutput
test_vectorOutput <- function(){
    #file <- system.file("extdata", "test4.csv", package="PhenStat")
    data_typical <- read.csv("PhenStat/inst/extdata/test4.csv")
    test <- PhenList(dataset=data_typical,testGenotype="Mysm1/+", dataset.clean=TRUE)
    result <- testDataset(test,depVariable="Lean.Mass")
    vector_results <- vectorOutput(result)   
    checkEquals(length(vector_results), 34, "length of vector is not 33") 
    checkTrue(!is.null(names(vector_results)),"names of vector are missed")
    checkTrue(all(lapply(vector_results,is.character)==TRUE),"not all vector elements are of type character")
}
###############################################################################
# generateGraphs for different stat. analysis methods (MM, TF, RR, FE)
# NB! Problem with 'cairo' argument (works only if cairo package is installed)
test_generateGraphs <- function(){
    #file <- system.file("extdata", "test4.csv", package="PhenStat")
    dirName <- ".graphics"
    # MM
    data_typical <- read.csv("PhenStat/inst/extdata/test4.csv")
    test <- PhenList(dataset=data_typical,testGenotype="Mysm1/+", dataset.clean=TRUE,outputMessages=FALSE)
    result <- testDataset(test,depVariable="Lean.Mass",outputMessages=FALSE)
    dir.create(dirName, showWarnings = FALSE)
    generateGraphs(result,dir=dirName,type="cairo")  
    checkEquals(length(list.files(dirName)), 7, "number of graphics for MM is not 7") 
    unlink(dirName, recursive = TRUE)
    ## TF
    file <- system.file("extdata", "test7_TFE.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="het",
            refGenotype = "WT",
            dataset.colname.sex="sex",
            dataset.colname.genotype="Genotype",
            dataset.values.female="f",
            dataset.values.male= "m",
            dataset.colname.weight="body.weight",
            dataset.colname.batch="Date_of_procedure_start",outputMessages=FALSE)
    test_TF <- TFDataset(test,depVariable="Cholesterol",outputMessages=FALSE)    
    # when "testDataset" function's argument "callAll" is set to FALSE 
    # only "startTFModel" function is called - the first step of TFE framework
    result_TF  <- testDataset(test_TF,
            depVariable="Cholesterol",
            method="TF",outputMessages=FALSE)
    dir.create(dirName, showWarnings = FALSE)
    generateGraphs(result_TF,dir=dirName,type="cairo")  
    checkEquals(length(list.files(dirName)), 5, "number of graphics for TF is not 5") 
    unlink(dirName, recursive = TRUE)
    ## RR
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
        testGenotype="Sparc/Sparc",outputMessages=FALSE)
    # "RRTest" function is called from "testDataset" function
    result_RR <- testDataset(test,
        depVariable="Lean.Mass",
        method="RR",outputMessages=FALSE)
    dir.create(dirName, showWarnings = FALSE)
    generateGraphs(result_RR,dir=dirName,type="cairo")  
    checkEquals(length(list.files(dirName)), 1, "number of graphics for RR is not 1") 
    unlink(dirName, recursive = TRUE)
    ## FE
    file <- system.file("extdata", "test_categorical.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
        testGenotype="Aff3/Aff3",outputMessages=FALSE)
    result_FE <- testDataset(test,depVariable="Thoracic.Processes",method="FE",outputMessages=FALSE)
    dir.create(dirName, showWarnings = FALSE)
    generateGraphs(result_FE,dir=dirName,type="cairo")  
    checkEquals(length(list.files(dirName)), 1, "number of graphics for RR is not 1") 
    unlink(dirName, recursive = TRUE)
}
###############################################################################
# test unit for TF function - lack of concurrent controls
test_TFOutput <- function(){
    data_typical <- read.csv("PhenStat/inst/extdata/test4.csv")
    test <- PhenList(dataset=data_typical,testGenotype="Mysm1/+", dataset.clean=TRUE,outputMessages=FALSE)
    checkException(testDataset(test,depVariable="Lean.Mass",method="TF"),
        "testDataset fails with check for sufficient data in TF framework")
    checkEquals(length(classificationTag(result)),1, "length of classificationTag output is not 1")
    checkTrue(is.character(classificationTag(result)), "classificationTag output is not character")
    checkEquals(classificationTag(result),
        "With phenotype threshold value 0.01 - significant in males and in combined dataset",
        "classificationTag output is incorrect for this dataset")
}
###############################################################################
# test unit for trasnformation function
test_transformation <- function(){
data1 <- read.csv("PhenStat/inst/extdata/permutationTest.csv")
test <- PhenList(dataset=data1,testGenotype="Kcne2/Kcne2", 
        dataset.colname.weight="Weight.Value",
        dataset.clean=TRUE,outputMessages=FALSE)
resultSodium <- testDataset(test,depVariable="Sodium",outputMessages=FALSE)
checkTrue(!(resultSodium@transformationRequired), "transformation is required")
resultSodiumAdded <- testDataset(test,depVariable="Sodium_added",outputMessages=FALSE)
checkTrue(resultSodiumAdded@transformationRequired, "transformation is not required")
checkEquals(resultSodiumAdded@lambdaValue,7.9, "transformation lambda value is not 7.9")
checkEquals(resultSodiumAdded@scaleShift,3, "transformation scale shift value is not 3")
resultHdl <- testDataset(test,depVariable="Hdl",outputMessages=FALSE)
checkEquals(resultHdl@lambdaValue,2.15, "transformation lambda value is not 2.15")
resultCa <- testDataset(test,depVariable="Ca",outputMessages=FALSE)
checkEquals(resultCa@lambdaValue,-2.9, "transformation lambda value is not -2.9")
resultTrigs<- testDataset(test,depVariable="Trigs",outputMessages=FALSE)
checkEquals(resultTrigs@lambdaValue,0, "transformation lambda value is not 0, which means log transformed")
#testDataset function with argument values requiring NOT TO TRANSFORM dependent variable
resultTrigs<- testDataset(test,depVariable="Trigs",transformValues=FALSE,outputMessages=FALSE)
checkTrue(!resultTrigs@transformationRequired, "transformation is not required")
resultAst <- testDataset(test,depVariable="Ast",outputMessages=FALSE)
checkEquals(resultAst@lambdaValue,0.2, "transformation lambda value is not 0.2")
}
###############################################################################

