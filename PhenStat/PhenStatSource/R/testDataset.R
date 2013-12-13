## Copyright © 2011-2013 EMBL - European Bioinformatics Institute
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
testDataset <- function(phenList, depVariable, equation="withWeight", 
        outputMessages=TRUE, pThreshold=0.05, method="MM", callAll=TRUE, keepList=NULL)
{
    
    stop_message <- ""
    
    ## Checks and stop messages       
    if(is(phenList,"PhenList")) {
        x <- phenList$dataset    
        
        if (!(depVariable %in% colnames(x)))
        stop_message <- paste("Error:\nDependent variable column '",
                depVariable,"' is missed in the dataset.\n",sep="")
        else {
            columnOfInterest <- x[,c(depVariable)]
            ## Test: depVariable is numeric in MM framework
            if(is.numeric(columnOfInterest)){
                variability <- length(unique(columnOfInterest))/length(columnOfInterest)
                ## Test: depVariable is continuous variable in MM framework (sufficcient variablity)   
                if ((variability<0.005) && method=="MM"){ 
                    stop_message <- paste(stop_message,"Error:\nInsufficient 
                            variability in the dependent variable '",depVariable,
                            "' for MM framework. Fisher Exact Test can be better way to 
                            do the analysis.\n",sep="") 
                }    
            }
            else if (method=="MM"){
                method <- "FE"
                if (outputMessages)
                message(paste("Warning:\nDependent variable '",depVariable,
                                "' is not numeric. Fisher Exact Test will be used for the 
                                analysis of this dependent variable.",sep=""))
            }
            
            if (method=="MM"){
                ## Test: depVariable variablity in Genotypes (require at least 2 levels)            
                Genotype_levels <- levels(x$Genotype)
                Gender_levels <- levels(x$Gender)  
                
                passed <- TRUE
                for (i in 1:length(Genotype_levels)){
                    GenotypeSubset <- subset(x, x$Genotype==Genotype_levels[i])
                    for (j in 1:length(Gender_levels)){           
                        GenotypeGenderSubset <- subset(GenotypeSubset, 
                                GenotypeSubset$Gender==Gender_levels[j]) 
                        
                        columnOfInterest <- GenotypeGenderSubset[,c(depVariable)]
                        if (length(unique(columnOfInterest))==1){
                            passed<-FALSE
                        }
                    }                
                }
                if (!passed)
                stop_message <- paste(stop_message,"Error:\nInsufficient 
                        variability in the dependent variable '",depVariable,
                        "' for genotype/gender combinations to allow the application 
                        of Mixed Model.\n",sep="")             
            }
            else if (method=="FE") {
                ## Test: depVariable number of levels is at least 2
                if (length(levels(factor(columnOfInterest,exclude=NA)))==0){
                    stop_message <- paste("Error:\nInsufficient data in the 
                            dependent variable '",depVariable,
                            "' to allow the application of Fisher Exact test framework.\n",sep="") 
                }
                if (length(levels(factor(columnOfInterest,exclude=NA)))>10){
                    stop_message <- paste("Error:\nPackage supports up to 10 
                            levels in dependent variable in FE framework. 
                            The variable '",depVariable,"' has more than 10 levels.\n",sep="") 
                }    
            }     
        }
        
    } else {
        stop_message <- "Error:\nPlease create a PhenList object first.\n"
    }    
    
    
    if (!(equation %in% c("withWeight","withoutWeight")) && method=="MM")
    stop_message <- paste(stop_message,"Error:\nPlease define equation you 
            would like to use from the following options: 'withWeight', 
            'withoutWeight'.\n",sep="")
    
    ## Dealing with provided significance values
    if (!is.null(keepList)){
        ## Stop function if there are no enough needed input parameters
        if ((length(keepList[keepList==TRUE]) + length(keepList[keepList==FALSE])) !=5) 
        stop_message <- paste(stop_message,
                "Error:\nPlease define the values for 'keepList' list, 
                where for each effect/part of the model TRUE/FALSE value defines to 
                keep it in the model or not, for instance: 
                'keepList=c(keepBatch=TRUE,keepEqualVariance=TRUE,keepWeight=FALSE,
                        keepGender=TRUE,keepInteraction=TRUE)'.\n",sep="")
        
    }
    
    if (nchar(stop_message)>0){
        if (outputMessages)   
        message(stop_message)
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))      
        stop()
    }
    
    if (!('Weight' %in% colnames(x)) && equation=="withWeight" && method=="MM"){
        if (outputMessages)
        message("Warning:\nWeight column is missed in the dataset. Equation 
                'withWeight' can't be used and has been replaced to 'withoutWeight'.\n")
        equation <- "withoutWeight"
    }
    
    if (depVariable=='Weight' && equation=="withWeight" && method=="MM"){
        if (outputMessages)
        message("Warning:\nWeight is used as dependent variable. Equation 
                'withWeight' can't be used in such a case and has been replaced to 'withoutWeight'.\n")
        equation <- "withoutWeight"
    }
    
    ## END Checks and stop messages
    
    if (outputMessages)
    message(paste("Information:\nDependent variable: '",
                    depVariable,"'.\n",sep="")) 
    
    ## Mixed Models framework
    if (method=="MM")    { 
        
        if (outputMessages)
        message(paste("Information:\nMethod: Mixed Model framework.\n",sep="")) 
        
        result <- startModel(phenList, depVariable, equation, 
                outputMessages, pThreshold, keepList)
        
        ## Perform all framework methods 
        if (callAll){
            if (outputMessages)
            message(paste("Information:\nPerform all MM framework stages: 
                            startModel and finalModel.\n",sep=""))
            result <- finalModel(result, outputMessages)
        }
        
    }
    else if (method=="FE") {
        ## Fisher Exact Test 
        if (outputMessages)
        message(paste("Information:\nMethod: Fisher Exact Test framework.\n",sep="")) 
        result <- FisherExactTest(phenList,depVariable,outputMessages)
    }
    else {
        if (outputMessages)   
        message(paste("Error:\nMethod define in the 'method' argument '",
                        method,"' is not supported.\nAt the moment we are supporting 'MM' 
                        value for Mixed Model framework and 'FE' value for Fisher Exact 
                        Test framework.\n",sep=""))
        
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))      
        stop()
        
    }        
    return(result)   
}
##------------------------------------------------------------------------------