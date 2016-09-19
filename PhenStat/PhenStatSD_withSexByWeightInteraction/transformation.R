## Copyright © 2012-2015 EMBL - European Bioinformatics Institute
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
## transformation.R contains functions for the Box-Cox transformation
##------------------------------------------------------------------------------
##Function transforms the given vector of values according to the lambda value:
# log transformation if the lambda is 0, power transformation otherwise
performTransformation <- function (values, lambda, scaleShift){

    
    if (lambda==0){
        transformedValues <- lapply(values, 
                function (x) {log(x+scaleShift)})
    }
    else {
        transformedValues <- lapply(values, 
                function (x) {y<-x+scaleShift;((sign(y)*abs(y)^lambda)-1)/lambda})
        #(((values+scaleShift)^lambda)-1)/lambda   
       
    }
    
    return(unlist(transformedValues))
}
##------------------------------------------------------------------------------
##Function reverse back the transformed values according to the lambda value:
# exponential transformation if the lambda is 0, fractional power transformation otherwise
performReverseTransformation <- function (values, lambda, scaleShift){
    if (lambda==0){
        reverseValues <- lapply(values, function (x) {exp(x) - scaleShift})
    }
    else {       
        reverseValues <- lapply(values, 
                function (x) {y<-((x * lambda)+1);(sign(y)*abs(y)^(1/lambda))-scaleShift})
        
    }
    return(unlist(reverseValues))
}
##------------------------------------------------------------------------------
##Function role
#1.  Rescale data if needed (Box-Cox transformation can only be applied on data >0)
#2.  Calculate lambda value for a Box-Cox transformation 
#3.  Assess whether transformation is required (if confidence interval includes 1 then no transformation required)
#4.  Return output to allow application of the transformation if necessary using 
#    the returned lambda value and any associated rescaling that was needed.
#5.  Requires a function that will be fitted determined from formulaAssessingBoxCox

#Argument form:
#multipleBatches Yes/No   
#noSexes  1/2
#depVariable string
#data - requires column Assay.Date, Gender, Genotype, Weight if weight included, 
#WeightIncluded - Yes/No

#Output form:
#Vector 5 elements - all numeric
#TransfomrationRequired TRUE,  FALSE
determiningLambda <- function(phenList, depVariable, equation="withWeight"){
    df <- getDataset(phenList)
    noSexes <- noSexes(phenList)
    multipleBatches <-  ifelse(multipleBatches(phenList), "Yes", "No")
    transformationCode <- 100
    if (!weightIn(phenList)){
        equation <- "withoutWeight"
    }
    #Box-Cox process requires all values to be greater than zero - we need to add a constant to shift the scale
    
    scaleShiftRequire=min(df[, depVariable], na.rm = TRUE)<=0
    
    if(scaleShiftRequire==TRUE){        
        #create a new column for scaled variable        
        newName=paste(depVariable, "scaled", sep="_")
        df[newName] <- NA
        
        #the rescaling has two elements - First you add the minimum value seen so that it is now zero 
        # and then you add 1 to move it away from zero
        df[newName]= df[ ,depVariable]+abs(min(df [,depVariable], na.rm = TRUE))+1    
        scaleShift=abs(min(df[ ,depVariable], na.rm = TRUE))+1    
        formulatoTest=formulaDeterminingLambda(noSexes, newName, multipleBatches, equation)
        
    }else{
        
        formulatoTest=formulaDeterminingLambda(noSexes, depVariable, multipleBatches, equation)
        scaleShift = NA    
    }
    
    #determine lambda using the formula returned
    boxcox_out=boxcox(formulatoTest, data = df, plotit = FALSE, lambda = seq(-100, 100, 1/10))
    #determine the 95% confidence interval for lambda
    lambda_CI=range(boxcox_out$x[boxcox_out$y > max(boxcox_out$y)-qchisq(0.95,1)/2])
    #determine the midepoint of the 95% confidence interval
    midpoint=round(boxcox_out$x[which.max(boxcox_out$y)],digits=4)  #(lambda_CI[2]+lambda_CI[1])/2
    
    #if the lambda value is close to zero then it reverts to a log transfomration.  
    # I have arbituarily chosen |0.1| as threshold (Natasha Karp). 
    # Based on observation that 0.5 equals a sqrt transformation and 0.2 did ok in a test dataset.    
    if (midpoint > -0.1 & midpoint < 0.1){
        lambda_value=0
        transformationCode <- 2
    }else{
        lambda_value=midpoint
        transformationCode <- 3
    }
    #determine if the confidence interval includes 1 where the transformation has no impact  
    #if includes one - TransformationRequired=FALSE if excludes one then TransformationRequired=TRUE
    #RangeExcludes1= !(1 > lambda_CI[1] & 1 < lambda_CI[2])
    RangeExcludes1 <- TRUE
    
    if (1 >= lambda_CI[1] & 1 <= lambda_CI[2]){
        transformationCode <- 1 
        RangeExcludes1 <- FALSE
    }
    
    if ((lambda_value > 5) || (lambda_value < -5)){
      transformationCode <- 4   
      RangeExcludes1 <- FALSE
    }

    #collate results for export
    output=c(lambda_CI, midpoint, lambda_value, as.character(RangeExcludes1), scaleShift, transformationCode)
    names(output)=c("CI_min", "CI_max", "midpointCI", "lambda_value", "TransformationRequired", 
            "scalingConstant", "Code")
    return(output)
}

#output is a vector with four elements.  The first three are numeric and the lsat is boolean with 1=tRUE and 2=FALSE

##------------------------------------------------------------------------------
##formulaDeterminingLambda returns the starting formula that you wish to interogate 
#the data with in estimating the lambda for a box cox transformation
#multipleBatches Yes/No   
#noSexes  1/2
#depVariable string
#WeightIncluded Yes/No
#returns: formula 
formulaDeterminingLambda <- function(noSexes, depVariable, multipleBatches, equation){
    
    formulaToTest <- switch(equation,    
            
            withWeight = {
                
                switch(multipleBatches, 
                        Yes= {
                            if(noSexes==2){
                                as.formula(paste(depVariable, "~", paste("Genotype","Sex", "Genotype*Sex", 
                                                        "Weight", "Batch", sep="+")))                            
                            }else{
                                as.formula(paste(depVariable, "~", paste("Genotype", "Weight", "Batch", sep="+")))
                            }
                        },
                        
                        No = {
                            if(noSexes==2){
                                as.formula(paste(depVariable, "~", paste("Genotype","Sex", "Genotype*Sex", 
                                                        "Weight", sep="+")))                            
                            }else{
                                as.formula(paste(depVariable, "~", paste("Genotype", "Weight", sep="+")))
                            }
                        }
                        )},
            withoutWeight = {
                
                switch(multipleBatches, 
                        Yes= {
                            if(noSexes==2){
                                as.formula(paste(depVariable, "~", paste("Genotype","Sex", 
                                                        "Genotype*Sex", "Batch", sep="+")))                            
                            }else{
                                as.formula(paste(depVariable, "~", paste("Genotype",  "Batch", sep="+")))
                            }
                        },
                        
                        No = {
                            if(noSexes==2){
                                as.formula(paste(depVariable, "~", paste("Genotype","Sex", "Genotype*Sex",  sep="+")))                            
                            }else{
                                as.formula(paste(depVariable, "Genotype", sep="~"))
                            }
                        })
            })
    return(formulaToTest)
}
##------------------------------------------------------------------------------