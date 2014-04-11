## Copyright � 2011-2013 EMBL - European Bioinformatics Institute
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
        outputMessages=TRUE, pThreshold=0.05, method="MM", callAll=TRUE, keepList=NULL, dataPointsThreshold=4)
{
    
    stop_message <- ""
    
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
    if (!(equation %in% c("withWeight","withoutWeight")) && method=="MM")
    stop_message <- paste(stop_message,"Error:\nPlease define equation you ", 
            "would like to use from the following options: 'withWeight', 'withoutWeight'.\n",sep="")
    
    # 4  Checks for rovided significance values
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
    if (!(method %in% c("MM","FE"))){
        stop_message <- paste(stop_message,"Error:\nMethod define in the 'method' argument '",
                method,"' is not supported.\nAt the moment we are supporting 'MM' ", 
                "value for Mixed Model framework and 'FE' value for Fisher Exact Test framework.\n",sep="")
        
    }   
    
    # 6
    if (dataPointsThreshold<2) {
        dataPointsThreshold <- 2
        if (outputMessages)
        message("Warning:\nData points threshold is set to 2 (minimal value).\n")
    }
    
    # 7 
    if (nchar(stop_message)==0) {
        x <- phenList$dataset 
        checkDepV <- columnChecks(x,depVariable,dataPointsThreshold)
        
        # Presence
        if (!checkDepV[1])
        stop_message <- paste("Error:\nDependent variable column '",
                depVariable,"' is missed in the dataset.\n",sep="")
    }
    
    
    ## STOP CHECK ARGUMENTS   
    
    ## DATASET'S CHECKS   
    # Dataset checks depending on selected method
    if (nchar(stop_message)==0){
        
        
        checkDepVLevels <- columnLevels(x,depVariable)   

        checkWeight <- columnChecks(x,"Weight",dataPointsThreshold)
        # MM checks
        # Should go first since there are switches between MM to FE
        if (method=="MM"){

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
                    stop_message <- paste(stop_message,"Error:\nNo variability in dependent variable '",
                            depVariable,"'.\n",sep="") 
                } 
                else  if (variability<0.005){ 
                    stop_message <- paste(stop_message,"Error:\nInsufficient variability in dependent variable '",
                            depVariable,
                            "' for MM framework. Fisher Exact Test can be better way to do the analysis.\n",sep="") 
                } 
                # Data points - number of data points in the depVariable column for genotype/sex combinations
                else if (!checkDepV[3])
                stop_message <- paste(stop_message,"Error:\nNot enough data points in dependent variable '",
                        depVariable,
                        "' for genotype/sex combinations to allow the application of Mixed Model. ",
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
                        
                        # Equality to depVariable
                        columnOfInterest <- na.omit(x[,c(depVariable)])
                        columnOfWeight <- na.omit(x$Weight)
                        
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
        
        result <- startModel(phenList, depVariable, equation, 
                outputMessages, pThreshold, keepList)
        
        ## Perform all framework methods 
        if (callAll && is(result,"PhenTestResult")){
            result <- finalModel(result, outputMessages)
        }
        
    }
    else if (method=="FE") {
        ## Fisher Exact Test 
        if (outputMessages)
        message(paste("Information:\nMethod: Fisher Exact Test framework.\n",sep="")) 
        result <- FisherExactTest(phenList,depVariable,outputMessages)
    }
    
    return(result)   
}
##------------------------------------------------------------------------------
## Returns values after null removing procedure: No of data points, No of levels, 
## No of Genotype/Sex combinations, No of data points for each one combination
columnLevels <- function(dataset, columnName){
    
    columnOfInterest <- na.omit(dataset[,c(columnName)])
    # Test if there are data points in a column other than NA
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