## Copyright © 2012-2014 EMBL - European Bioinformatics Institute
## 
## Licensed under the Apache License, Version 2.0 (the "License"); 
## you may not use this file except in compliance with the License.  
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##------------------------------------------------------------------------------
## testDataset.R contains testDataset and modelFormula functions and constructs 
## the PhenTestResult class object
##------------------------------------------------------------------------------
## testDataset function performs the following checks for dependent variable:
## 1. "depVariable" column should present in the dataset
## 2. "depVariable" should be numeric for Mixed Model (MM) framework, 
## otherwise performs Fisher Exact Test (FE)
## 3. Each one genotype level should have more than one "depVariable" level
## (variability) for MM framework, otherwise recommends FE framework. 
## MM framework consists of two functions: startModel and finalModel. 
## callAll argument set to TRUE (default) indicates to call both functions
## FET framework consists of only one function FisherExactTest that creates 
## count matrix (matrices) and perform test(s).  
##------------------------------------------------------------------------------
testDataset <- function(phenList=NULL, depVariable=NULL, equation="withWeight", 
        outputMessages=TRUE, pThreshold=0.05, method="MM", callAll=TRUE, 
        keepList=NULL, dataPointsThreshold=4, RR_naturalVariation=95, 
        RR_controlPointsThreshold=60, transformValues=TRUE)
{

    stop_message <- ""
    transformationRequired <- FALSE
    lambdaValue <- NA
    columnOfBatch <- NULL
    columnOfWeight <- NULL
    columnOfInterestAdjusted <- NULL
    ## CHECK ARGUMENTS   
    
    # 1
    if (is.null(phenList) || !is(phenList,"PhenList")){
        stop_message <- paste(stop_message,"Error:\nPlease create and specify PhenList object.\n",sep="")
    }
    
    # 2
    if (is.null(depVariable)){
        stop_message <- paste(stop_message,"Error:\nPlease specify dependent variable, ",
                "for example: depVariable='Lean.Mass'.\n",sep="")
    }
    
    # 3
    if (!(equation %in% c("withWeight","withoutWeight")) && method %in% c("MM","TF"))
        stop_message <- paste(stop_message,"Error:\nPlease define equation you ", 
            "would like to use from the following options: 'withWeight', 'withoutWeight'.\n",sep="")
    
    # 4  Checks for provided significance values
    if (!is.null(keepList)){
        ## Stop function if there are no enough needed input parameters
        if ((length(keepList[keepList==TRUE]) + length(keepList[keepList==FALSE])) !=5) 
        stop_message <- paste(stop_message,
                "Error:\nPlease define the values for 'keepList' list ", 
                "where for each effect/part of the model TRUE/FALSE value defines to ",
                "keep it in the model or not, for instance: ", 
                "'keepList=c(keepBatch=TRUE,keepEqualVariance=TRUE,keepWeight=FALSE, ",
                "keepSex=TRUE,keepInteraction=TRUE)'.\n",sep="")
        
    }
    
    # 5
    if (!(method %in% c("MM","FE","RR","TF","LR"))){
        stop_message <- paste(stop_message,"Error:\nMethod define in the 'method' argument '",
                method,"' is not supported.\nAt the moment we are supporting 'MM' ", 
                "value for Mixed Model framework, 'FE' value for Fisher Exact Test framework",
                "'RR' for Reference Ranges Plus framework and 'TF' value for Time as Fixed Effect framework.\n",
                sep="")
        
    }   
    
    # 6
    if (dataPointsThreshold<2) {
        dataPointsThreshold <- 2
        if (outputMessages)
        message("Warning:\nData points threshold is set to 2 (minimal value).\n")
    }
    
    # RR_naturalVariation=95, RR_controlPointsThreshold=60
    if (RR_naturalVariation<60) {
        RR_naturalVariation <- 60
        if (outputMessages)
        message("Warning:\nNatural variation threshold is set to 60 (minimal value).\n")
    }
    
    if (RR_controlPointsThreshold<40) {
        RR_controlPointsThreshold <- 40
        if (outputMessages)
        message("Warning:\nControl points threshold is set to 40 (minimal value).\n")
    }
  
    # 7 
    if (nchar(stop_message)==0) {
        # create small data frame with values to analyse
        columnOfInterest <- getColumn(phenList,depVariable)
        # Presence
        if (is.null(columnOfInterest))
            stop_message <- paste("Error:\nDependent variable column '",
                depVariable,"' is missed in the dataset.\n",sep="")
    }    
    
    ## STOP CHECK ARGUMENTS   
    
    # Creates a subset with analysable variable
    if (nchar(stop_message)==0) {
            checkDepV <- columnChecks(getDataset(phenList),depVariable,dataPointsThreshold)
            # TRANSFORMATION
            # check if the method selected is not for categorical data - if so we don't want to transform values
            if (!(method %in% c("FE","LR","RR")) && checkDepV[2] && checkDepV[3]) {
                # check for transformation
                transformationVector <- determiningLambda(phenList,depVariable,equation)
                transformationRequired <- as.logical(transformationVector[[5]])
                lambdaValue <- as.numeric(transformationVector[[4]])
                if (!is.na(transformationVector[[6]])){
                    scaleShift <- as.numeric(transformationVector[[6]])
                }
                else {scaleShift <- 0 }
                if (transformValues && transformationRequired){
                    columnOfInterestOriginal <- columnOfInterest
                    columnOfInterest <- round(transformValues(columnOfInterest,lambdaValue,scaleShift),digits=6)
                }
                if (batchIn(phenList)){
                    # Adjusted for batch depVariable values WITHOUT transformation!
                    columnOfInterestAdjusted=getColumnBatchAdjusted(phenList,depVariable)
                }
            }
            else {
                columnOfInterestOriginal <- columnOfInterest
                # Transformation for LR -> 0/1    
            }


            columnOfSex <- getColumn(phenList,"Sex")
            columnOfGenotype <- getColumn(phenList,"Genotype")
            if (batchIn(phenList)){
                columnOfBatch <- getColumn(phenList,"Batch")
            }

            if (weightIn(phenList)){
                columnOfWeight <- getColumn(phenList,"Weight")
            }
            # Get rid of factors 
            if (class(columnOfInterest)=="factor"){
                columnOfInterest <- as.character(columnOfInterest)
                tryCatch({
                            # try to convert into numbers
                            columnOfInterest<-as.numeric(columnOfInterest)
                        },
                        warning = function(war){
                            # convert into characters
                            columnOfInterest <- as.character(columnOfInterest)
                        },
                        error = function(err){
                            # convert into characters
                            columnOfInterest <- as.character(columnOfInterest)        
                        }) 
            }
            
            if (!is.null(columnOfBatch) && !is.null(columnOfWeight)){
                datasetToAnalyse <- data.frame(columnOfInterest,columnOfSex,columnOfGenotype,columnOfBatch,columnOfWeight)
                colnames(datasetToAnalyse)<-c(depVariable,"Sex","Genotype","Batch","Weight")       
            }
            else if (!is.null(columnOfBatch)){
                datasetToAnalyse <- data.frame(columnOfInterest,columnOfSex,columnOfGenotype,columnOfBatch)
                colnames(datasetToAnalyse)<-c(depVariable,"Sex","Genotype","Batch")    
            } 
            else if (!is.null(columnOfWeight)){
                    datasetToAnalyse <- data.frame(columnOfInterest,columnOfSex,columnOfGenotype,columnOfWeight)
                    colnames(datasetToAnalyse)<-c(depVariable,"Sex","Genotype","Weight")  
                }
            else {
                datasetToAnalyse <- data.frame(columnOfInterest,columnOfSex,columnOfGenotype)
                colnames(datasetToAnalyse)<-c(depVariable,"Sex","Genotype") 
                }
            if (transformValues && transformationRequired && (!(method %in% c("FE","LR","RR")))){
                columnNameOriginal <- paste(depVariable,"_original",sep="")
                datasetToAnalyse[,columnNameOriginal] <- columnOfInterestOriginal
            }
            if (!is.null(columnOfInterestAdjusted)) {
                columnNameAdjusted <- paste(depVariable,"_adjusted",sep="")
                datasetToAnalyse[,columnNameAdjusted] <- columnOfInterestAdjusted  
            }


        phenListToAnalyse <- new("PhenList",datasetPL=datasetToAnalyse,
                refGenotype = refGenotype(phenList),
                testGenotype = testGenotype(phenList),
                hemiGenotype = hemiGenotype(phenList))
        
        x <- datasetToAnalyse
        
        checkDepV <- columnChecks(x,depVariable,dataPointsThreshold)
    
        checkDepVLevels <- columnLevels(x,depVariable)   

        checkWeight <- columnChecks(x,"Weight",dataPointsThreshold)

        # MM checks
        # Go first since there are switches between MM to FE
        if (method %in% c("MM","TF")){

            # Numeric
            if (!checkDepV[2]) {
                method <- "FE"
                if (outputMessages)
                message(paste("Warning:\nDependent variable '",depVariable,
                                "' is not numeric. Fisher Exact Test will be used for the ", 
                                "analysis of this dependent variable.\n",sep=""))
            }
            else{
                
               
                
                # Variability - the ratio of different values to all values in the column
                if (checkDepVLevels[1]>0)
                    variability <- checkDepVLevels[2]/checkDepVLevels[1] 
                else 
                    variability <- 0
                # where checkDepVLevels[2] contains number of levels and checkDepVLevels[1] contains number of data points
                # One level only
                if (checkDepVLevels[2]==1 || checkDepVLevels[2]==0){ 
                    if (transformValues && transformationRequired && (!(method %in% c("FE","LR","RR")))){
                         stop_message <- paste(stop_message,"Error:\nNo variability in dependent variable '",
                                        depVariable,". Try with argument transformValues set to FALSE.'.\n",sep="") 
                    }
                    else {    
                         stop_message <- paste(stop_message,"Error:\nNo variability in dependent variable '",
                                        depVariable,"'.\n",sep="") 
                    }
                } 
                else  if (variability<0.005){ 
                    stop_message <- paste(stop_message,"Error:\nInsufficient variability in dependent variable '",
                            depVariable,
                            "' for MM/TF framework. Fisher Exact Test can be better way to do the analysis.\n",sep="") 
                } 
                # Data points - number of data points in the depVariable column for genotype/sex combinations
                else if (!checkDepV[3])
                stop_message <- paste(stop_message,"Error:\nNot enough data points in dependent variable '",
                        depVariable,
                        "' for genotype/sex combinations to allow the application of MM/TF framework. ",
                        "Threshold used: ",dataPointsThreshold,".\n",sep="")     
                
                #Weight checks
                if (equation=="withWeight") {
                    if (! checkWeight[1]){
                        if (outputMessages)
                        message(paste("Warning:\nWeight column is not present in dataset.\n",
                                        "Equation 'withWeight' can't be used in such a case and has been ",
                                        "replaced to 'withoutWeight'.\n", sep=""))  
                        equation <- "withoutWeight"   
                        
                    }
                    else{
                        #CHECKS FOR WEIGHT
                        # Equality to depVariable
                        columnOfInterest <- na.omit(columnOfInterest)
                        columnOfWeight <- na.omit(columnOfWeight)
                        
                        if (length(columnOfInterest) == length(columnOfWeight))
                            
                        if (sum(columnOfInterest-columnOfWeight) == 0){    
                            if (outputMessages)
                            message(paste("Warning:\nWeight and dependent variable values seemed to be equivalent. ",
                                            "Equation 'withWeight' can't be used in such a case and ",
                                            "has been replaced to 'withoutWeight'.\n", sep=""))
                            equation <- "withoutWeight"
                        }
                        
                        if (! checkWeight[2]){
                            if (outputMessages)
                            message(paste("Warning:\nWeight column is not numeric.\n",
                                            "Equation 'withWeight' can't be used in such a case and has been ",
                                            "replaced to 'withoutWeight'.\n", sep=""))
                            equation <- "withoutWeight"       
                        }    
                        
                        if (! checkWeight[3]){
                            if (outputMessages)
                            message(paste("Warning:\nWeight column does not have enough data points ",
                                            "for for genotype/sex combinations.\n",
                                            "Equation 'withWeight' can't be used in such a case and has been ",
                                            "replaced to 'withoutWeight'.\n", sep="")) 
                            equation <- "withoutWeight"   
                        } 
                                
                    }     
                }
            }            
        }
        
        # TF checks
        if (method=="TF" && nchar(stop_message)==0){
            # Lists possible combinations
            Genotype_levels <- levels(factor(x$Genotype))
            Batch_levels <- levels(factor(x$Batch))
            if (length(Batch_levels)<2 || length(Batch_levels)>5){
                stop_message <- paste(stop_message,"Error:\n'TF' framework requires from 2 to 5 batch levels. ",
                        "There is/are '",
                        length(Batch_levels),"' batch level(s) in the dataset.\n",sep="") 
            }
            TFDataPoints <- TRUE
            for (i in 1:length(Batch_levels)){
                BatchSubset <- subset(x, x$Batch==Batch_levels[i])
                for (j in 1:length(Genotype_levels)){           
                    GenotypeBatchSubset <- subset(BatchSubset, 
                            BatchSubset$Genotype==Genotype_levels[j]) 
                    columnOfInterestSubset <- na.omit(GenotypeBatchSubset[,c(depVariable)])                    
                    if (length(columnOfInterestSubset)==0){
                         TFDataPoints <- FALSE
                    }
                }   
   
            }
            if (!TFDataPoints)
                stop_message <- paste(stop_message,"Error:\n'TF' framework requires data points for at least one sex ",
                    "in all genotype/batch level combinations.\n",sep="") 
        }
        
        # FE checks
        if (method=="FE"){
            if (checkDepVLevels[2]==0)
            stop_message <- paste("Error:\nInsufficient data in the dependent variable '",
                    depVariable,
                    "' to allow the application of Fisher Exact test framework.\n",sep="") 
            if (checkDepVLevels[2]>10)
            stop_message <- paste("Error:\nPackage supports up to 10 levels ",
                    "in dependent variable in FE framework. ", 
                    "The variable '",depVariable,"' has more than 10 levels.\n",sep="") 
        } 

        # LR checks
        if (method=="LR"){
            if (checkDepVLevels[2]==0)
                stop_message <- paste("Error:\nInsufficient data in the dependent variable '",
                    depVariable,
                    "' to allow the application of Logistic Regression framework.\n",sep="") 
            if (checkDepVLevels[2]!=2){
                stop_message <- paste("Error:\nLogistic regression is applicable only if there are two levels ",
                    "in dependent variable. ", 
                    "The variable '",depVariable,"' has less or more than 2 levels.\n",sep="") 
            }
            else {
                depVariableLevels <- levels(factor(na.omit(columnOfInterest)))
                if (!(depVariableLevels[1]==0 && depVariableLevels[2]==1))
                    stop_message <- paste("Error:\nLogistic regression is applicable only if there are two levels ",
                    "in dependent variable: 0/1\n",sep="")                
            }
        } 
        
        # RR checks
        if (method=="RR"){
            # Numeric
            if (!checkDepV[2]) {
                method <- "FE"
                if (outputMessages)
                message(paste("Warning:\nDependent variable '",depVariable,
                                "' is not numeric. Fisher Exact Test will be used for the ", 
                                "analysis of this dependent variable.\n",sep=""))
            }
            else{
            
                if (checkDepVLevels[2]==0)
                    stop_message <- paste("Error:\nInsufficient data in the dependent variable '",
                        depVariable,
                        "' to allow the application of RR plus framework.\n",sep="") 
                # Number of control data
                
                
                controlSubset <- subset(x, x$Genotype==refGenotype(phenList))
                columnOfInterestSubset <- na.omit(controlSubset[,c(depVariable)])
                
                Sex_levels <- levels(factor(x$Sex))
                
                controlNotEnough <- FALSE
                errorText <- ""
                
                for (j in 1:length(Sex_levels)){           
                    GenotypeSexSubset <- subset(controlSubset, 
                            controlSubset$Sex==Sex_levels[j]) 
                        
                    columnOfInterestSubset <- na.omit(GenotypeSexSubset[,c(depVariable)])
                    
                    errorText <- paste(errorText, Sex_levels[j], " - ", length(columnOfInterestSubset),", ",sep="") 
                    
                    if (length(columnOfInterestSubset)<RR_controlPointsThreshold) {
                        controlNotEnough <- TRUE
                    }
                }                
                errorText <- substr(errorText, 1, nchar(errorText)-2) 
                
                
                if (controlNotEnough) 
                    stop_message <- paste("Error:\nInsufficient data in the dependent variable '",
                        depVariable,
                        "' control subset (",errorText,
                                ") to allow the application of RR plus framework.",
                        "\nThe threshold is ",RR_controlPointsThreshold," datapoints. \n",sep="") 
            }
        } 
        
    }    
    
    ## STOP DATASET'S CHECKS 
    
    
    # If problems are deteckted
    if (nchar(stop_message)>0){
        if (outputMessages)  { 
            message(stop_message)
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt))      
            stop()
        }
        else {
            stop(stop_message)
        }
    }
    
    
    
    # RUN FRAMEWORK
    if (outputMessages)
    message(paste("Information:\nDependent variable: '",
                    depVariable,"'.\n",sep="")) 
    
    ## Mixed Models framework
    if (method=="MM")    { 
        if (callAll){
            if (outputMessages)
            message(paste("Information:\nPerform all MM framework stages: startModel and finalModel.\n",sep=""))
        }    
        
        
        
        if (outputMessages)
        message(paste("Information:\nMethod: Mixed Model framework.\n",sep="")) 
        
        result <- startModel(phenListToAnalyse, depVariable, equation, 
                outputMessages, pThreshold, keepList)
        
        ## Perform all framework methods 
        if (callAll && is(result,"PhenTestResult")){
            result <- finalModel(result, outputMessages)
        }
        
    }
    else if (method=="TF") {
        if (callAll){
            if (outputMessages)
            message(paste("Information:\nPerform all TF framework stages: startTFModel and finalTFModel.\n",sep=""))
        }    
        
        
        
        if (outputMessages)
        message(paste("Information:\nMethod: Time as Fixed Effect framework.\n",sep="")) 
        
        result <- startTFModel(phenListToAnalyse, depVariable, equation, 
                outputMessages, pThreshold, keepList)
        
        ## Perform all framework methods 
        if (callAll && is(result,"PhenTestResult")){
            result <- finalTFModel(result, outputMessages)
        }
    }
    else if (method=="FE") {
        ## Fisher Exact Test 
        if (outputMessages)
        message(paste("Information:\nMethod: Fisher Exact Test framework.\n",sep="")) 
        result <- FisherExactTest(phenListToAnalyse,depVariable,outputMessages)
    }
    
    else if (method=="RR"){
        ## RR Plus
        if (outputMessages)
        message(paste("Information:\nMethod: Reference Ranges Plus framework.\n",sep="")) 
        result <- RRTest(phenListToAnalyse,depVariable,outputMessages,RR_naturalVariation,RR_controlPointsThreshold)
    }
    else if (method=="LR"){
        ## Logistic Regression
        if (callAll){
            if (outputMessages)
                message(paste("Information:\nPerform all LR framework stages: startLRModel and finalLRModel.\n",sep=""))
        }    
            
        if (outputMessages)
        message(paste("Information:\nMethod: Logistic Regression framework.\n",sep="")) 
            
        result <- startLRModel(phenListToAnalyse, depVariable, outputMessages=TRUE, pThreshold=0.05)
        
        ## Perform all framework methods 
        if (callAll && is(result,"PhenTestResult")){
                result <- finalLRModel(result, outputMessages)
        }
    }

    if (transformValues && transformationRequired && (method!="FE") && (method!="LR")){
        result@transformationRequired <- as.logical(transformationRequired)
        result@lambdaValue <- lambdaValue
        result@scaleShift <- scaleShift
    }
    else {
        result@transformationRequired <- FALSE
    } 
    
    return(result)   
}
##------------------------------------------------------------------------------
## Returns values after null removing procedure: No of data points, No of levels, 
## No of Genotype/Sex combinations, No of data points for each one combination
columnLevels <- function(dataset, columnName){
    
    columnOfInterest <- na.omit(dataset[,c(columnName)])

        
    values<- c(length(columnOfInterest))
    
    #Test for the data points quantity for Genotype/sex combinations
    Genotype_levels <- levels(factor(dataset$Genotype))
    Sex_levels <- levels(factor(dataset$Sex))
    values<-append(values,length(levels(factor(columnOfInterest))))
    
    values<-append(values,length(Genotype_levels)*length(Sex_levels))
    
    for (i in 1:length(Genotype_levels)){
        GenotypeSubset <- subset(dataset, dataset$Genotype==Genotype_levels[i])
        for (j in 1:length(Sex_levels)){           
            GenotypeSexSubset <- subset(GenotypeSubset, 
                    GenotypeSubset$Sex==Sex_levels[j]) 
            
            columnOfInterestSubset <- na.omit(GenotypeSexSubset[,c(columnName)])

            values<-append(values,length(columnOfInterestSubset))
            
        }                
    }
    return (values)
}
##------------------------------------------------------------------------------
## Checks the column for eligibility, returns values: 
## presence of column, all data are numeric, No of levels passed checks (number of data points for each genotype/sex
## combination is at least equals to threshold)
columnChecks <- function(dataset, columnName, dataPointsThreshold=4){
    presence <- TRUE
    numeric <- FALSE
    levelsCheck <- 0
    variabilityThreshold <- 10
    # Test: dependent variable presence 
    if (!(columnName %in% colnames(dataset))){
        presence <- FALSE
    }    
    else {
        columnOfInterest <- na.omit(dataset[,c(columnName)])

        if(all(sapply(columnOfInterest,is.numeric))){
            numeric <- TRUE
        }

        dataPointsSummary <- columnLevels(dataset,columnName)
        
        NoCombinations <- dataPointsSummary[3]
        #message(dataPointsSummary[3])
        variabilityThreshold <- NoCombinations
        #if (NoCombinations==4)
        #variabilityThreshold <- 3 
        
        for (i in 1:NoCombinations){
            if (dataPointsSummary[3+i] >= dataPointsThreshold)
            levelsCheck <- levelsCheck+1
            
        }
        
    }
    
    values <- c(presence, numeric, (levelsCheck>=variabilityThreshold)) 
    
    return (values)
}

##------------------------------------------------------------------------------
recommendMethod <- function(phenList=NULL, depVariable=NULL, 
        outputMessages=TRUE)
{
    
# evaluate    dataPointsThreshold, RR_controlPointsThreshold
    stop_message <- ""
    dataPointsThreshold <- 4
    RR_controlPointsThreshold <- 40
    suggestedFramework <- "NO"
    
    ## CHECK ARGUMENTS   
    
    # 1
    if (is.null(phenList) || !is(phenList,"PhenList")){
        stop_message <- paste(stop_message,"Error:\nPlease create and specify PhenList object.\n",sep="")
    }
    
    # 2
    if (is.null(depVariable)){
        stop_message <- paste(stop_message,"Error:\nPlease specify dependent variable, ",
                "for example: depVariable='Lean.Mass'.\n",sep="")
    } 
    else {
        columnOfInterest <- getColumn(phenList,depVariable)
    }
    
    # stop if there is something wrong with the arguments
    if (nchar(stop_message)==0) {
        x <- getDataset(phenList)
        checkDepV <- columnChecks(x,depVariable,dataPointsThreshold)
        
        # Presence
        if (!checkDepV[1])
        stop_message <- paste("Error:\nDependent variable column '",
                depVariable,"' is missed in the dataset.\n",sep="")
    }
    
    
    # If problems are deteckted
    if (nchar(stop_message)>0){
        if (outputMessages)  { 
            message(stop_message)
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt))      
            stop()
        }
        else {
            stop(stop_message)
        }
    }
    
    
    ## STOP CHECK ARGUMENTS     
    

                
        checkDepVLevels <- columnLevels(x,depVariable) 
        


        # check for FE        
        if (checkDepVLevels[2]>0 && checkDepVLevels[2]<=10) {
            suggestedFramework <- "FE"
        }
        
        # check for LR        
        if (checkDepVLevels[2]==2) {
            suggestedFramework <- paste(suggestedFramework,", LR",sep="") 
        }
        
        
        if (checkDepV[2]) { # NUMERIC
        
            # VARIABILITY
            variabilityPass <- TRUE
        
            # Variability - the ratio of different values to all values in the column
            if (checkDepVLevels[1]>0)
                variability <- checkDepVLevels[2]/checkDepVLevels[1] 
            else 
                variability <- 0
            # where checkDepVLevels[2] contains number of levels and checkDepVLevels[1] contains number of data points
            # One level only
            if (checkDepVLevels[2]==1 || checkDepVLevels[2]==0){ 
                variabilityPass <- FALSE
            } 
            else  
                if (variability<0.005){ 
                    variabilityPass <- FALSE
                } 
                # Data points - number of data points in the depVariable column for genotype/sex combinations
                else if (!checkDepV[3])
                    variabilityPass <- FALSE
            # VARIABILITY
         
         
            # check for TF
            if ('Batch' %in% colnames(x) && variabilityPass){
                phenListTF <- TFDataset(phenList,depVariable,outputMessages=FALSE,forDecisionTree=FALSE)
                xTF <- getDataset(phenListTF)            
                
                # check for batches - shoud be from 2 to 5 batches
                if (length(levels(factor(xTF$Batch))) >= 2 && length(levels(factor(xTF$Batch))) <= 5) {   
                    TF <- TRUE 
                    Genotype_levels <- levels(factor(x$Genotype))
                    Batch_levels <- levels(factor(x$Batch))
                    # check for data points in all genotype/batch level combinations (records at least in on Sex)
                    for (i in 1:length(Batch_levels)){
                        BatchSubset <- subset(x, x$Batch==Batch_levels[i])
                        # Genotype loop
                        for (j in 1:length(Genotype_levels)){           
                            GenotypeBatchSubset <- subset(BatchSubset, 
                                    BatchSubset$Genotype==Genotype_levels[j]) 
                            columnOfInterestSubset <- na.omit(GenotypeBatchSubset[,c(depVariable)])
                            if (length(columnOfInterestSubset)==0){
                                TF <- FALSE
                            }
                        }                        
                    }
                    if (suggestedFramework=="NO"){
                        suggestedFramework <- "TF"
                    }
                    else {
                        suggestedFramework <- paste(suggestedFramework,", TF",sep="") 
                    }
                }
                
            }
            
            # check for MM
            if (variabilityPass) {     
                if (suggestedFramework=="NO"){
                    suggestedFramework <- "MM"
                }
                else {
                    suggestedFramework <- paste(suggestedFramework,", MM",sep="") 
                }                
            }
                
            # checks for RR
            controlNotEnough <- FALSE
            if (checkDepVLevels[2]!=0) {               
                controlSubset <- subset(x, x$Genotype==refGenotype(phenList))
                columnOfInterestSubset <- na.omit(controlSubset[,c(depVariable)])                   
                Sex_levels <- levels(factor(x$Sex))                                  
                for (j in 1:length(Sex_levels)){           
                    GenotypeSexSubset <- subset(controlSubset, 
                    controlSubset$Sex==Sex_levels[j])                        
                    columnOfInterestSubset <- na.omit(GenotypeSexSubset[,c(depVariable)])                       
                    if (length(columnOfInterestSubset)<RR_controlPointsThreshold) {
                        controlNotEnough <- TRUE
                    }
                }  
            } 
            else {
                controlNotEnough <- TRUE
            }   
                   
            if (!controlNotEnough) { 
                if (suggestedFramework=="NO"){
                    suggestedFramework <- "RR"
                }
                else {
                    suggestedFramework <- paste(suggestedFramework," and RR",sep="") 
                }
            } 
          
            }
        
        return(suggestedFramework)   
       
}
##------------------------------------------------------------------------------                 "