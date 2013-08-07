# Copyright © 2011-2013 EMBL - European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License"); 
# you may not use this file except in compliance with the License.  
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#-----------------------------------------------------------------------------------
# testDataset.R contains testDataset and modelFormula functions and constructs the PhenTestResult 
# class object
#-----------------------------------------------------------------------------------
testDataset <- function(phenList, depVariable, equation="withWeight", outputMessages=TRUE, pThreshold=0.05, method="MM", callAll=TRUE, keepList=NULL)
# Performs the following checks for dependent variable:
# 1. "depVariable" column should present in the dataset
# 2. "depVariable" should be numeric for Mixed Model (MM) framework, otherwise performs Fisher Exact Test (FE)
# 3. Each one genotype level should have more than one "depVariable" level
# (variability) for MM framework, otherwise recommends FE framework. 


# MM framework consists of two functions: startModel and finalModel. callAll argument set to TRUE (default) indicates
# to call both functions

# FET framework consists of only one function FisherExactTest that creates ount matrix (matrices) and perform test(s).  

{
    
    stop_message <- ""
    
    # Checks and stop messages
    
    if(is(phenList,"PhenList")) {
        x <- phenList$dataset    
        
    } else {
        stop_message <- "Error:\nPlease create a PhenList object first.\n"
    }    
    
    if (!(depVariable %in% colnames(x)))
    stop_message <- paste("Error:\nDependent variable column '",depVariable,"' is missed in the dataset.\n",sep="")
    
    
    if (!(equation %in% c("withWeight","withoutWeight")) && method=="MM")
    stop_message <- "Error:\nPlease define equation you would like to use from the following options: 'withWeight', 'withoutWeight'\n."
    
    if (!('Weight' %in% colnames(x)) && equation=="withWeight" && method=="MM"){
        if (outputMessages)
        message("Warning:\nWeight column is missed in the dataset. Equation 'withWeight' can't be used and has been replaced to 'withoutWeight'.")
        equation="withoutWeight"
    }
    
    # Test: depVariable is continuous variable
    
    columnOfInterest <- x[,c(depVariable)]
    
    if(is.numeric(columnOfInterest)){
        if ((length(unique(columnOfInterest))/length(columnOfInterest)<0.05) && outputMessages && method=="MM") 
        message(paste("Warning: Dependent variable '",depVariable,"' is numeric but seemed to be categorical because there is little variation. Fisher Exact Test can be better way to do the analysis than Mixed Models.\n",sep="")) 
    }
    else if (method=="MM"){
        method="FE"
        if (outputMessages)
        message(paste("Warning:\nDependent variable '",depVariable,"' is not numeric. Fisher Exact Test will be used for the analysis of this dependent variable.\n",sep=""))
    }
    # Test: depVariable variablity in Genotypes (require at least 2 levels)
    
    Genotype_levels=levels(x$Genotype)
    Gender_levels=levels(x$Gender)  
    
    for (i in 1:length(Genotype_levels)){
        GenotypeSubset <- subset(x, x$Genotype==Genotype_levels[i])
        for (j in 1:length(Gender_levels)){           
            GenotypeGenderSubset <- subset(GenotypeSubset, GenotypeSubset$Gender==Gender_levels[j]) 
            columnOfInterest <- GenotypeGenderSubset[,c(depVariable)]
            if (length(unique(columnOfInterest))==1)
            stop_message <- paste("Error:\nInsufficient variability in the dependent variable '",depVariable,"' for genotype/gender combinations to allow the application of Mixed Model or Fisher Exact test framework.\n",sep="")
        }                
    }  
    
    # Dealing with provided significance values
    if (!is.null(keepList)){
        # Stop function if there are no enough needed input parameters
        if (length(keepList)!=5) 
        stop_message <- "Error:\nPlease define the values for 'keepList' list, where for each effect/part of the model TRUE/FALSE value defines to keep it in the model or not: 
        'keepList=c(keepBatch,keepVariance,keepWeight,keepGender,keepInteraction)'.\n"
        
    }
    
    
    if (nchar(stop_message)>1){
        if (outputMessages)   
        message(stop_message)
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))      
        stop()
    }
    
    # END Checks and stop messages
    
    if (outputMessages)
        message(paste("Information:\nDependent variable: '",depVariable,"'.\n",sep="")) 
    
    # Mixed Models framework
    if (method=="MM")    { 
        
        if (outputMessages)
        message(paste("Information:\nMethod: Mixed Model framework.\n",sep="")) 
        
        result <- startModel(phenList, depVariable, equation, outputMessages, pThreshold, keepList)
        
        # Perform all framework methods 
        if (callAll){
            if (outputMessages)
            message(paste("Information:\nPerform all MM framework methods: startModel and finalModel\n",sep=""))
            result <- finalModel(phenList, result, outputMessages)
        }
        
    }
    else if (method=="FE") {
        # Fisher Exact Test 
        if (outputMessages)
            message(paste("Information:\nMethod: Fisher Exact Test framework.\n",sep="")) 
        
        result <- FisherExactTest(phenList,depVariable,outputMessages)
    }
    else {
        if (outputMessages)   
            message(paste("Error:\nMethod define in the 'method' argument '",method,"' is not supported.\nAt the moment we are supporting 'MM' value for Mixed Model framework and 'FE' value for Fisher Exact Test framework.\n",sep=""))
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))      
        stop()
        
    }        
    
    return(result)   
}

  

