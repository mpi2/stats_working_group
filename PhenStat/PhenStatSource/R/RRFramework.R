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
## RRFramework.R contains ... functions for refernce range plus analysis 
##------------------------------------------------------------------------------
# Test for normal distirbution
# Gaussian distribution - the reference interval is calculated as th emean +-2 SD
# Non-Gaussian distribution - percentiles top 97.5 and 
# bottom 2.5 percentiles are used to limit of the reference range. 
# sufficient numbers of samples are required (minimum 50, ideal 120)

RRTest <- function(phenList, depVariable, 
        outputMessages=TRUE, rangeThreshold=60, naturalVariation=95)
{
    x <- phenList$dataset
    
    RR_left_male <- NA
    RR_right_male <- NA
    valuesInRange_male <- NA
    allValues_male <- NA
    percentageOut_male <- NA
    RR_left_female <- NA
    RR_right_female <- NA
    valuesInRange_female <- NA
    allValues_female <- NA
    percentageOut_female <- NA

    numberofgenders <- length(levels(x$Gender))
    
    # Defined control subset
    GenotypeSubset <- subset(x, x$Genotype==phenList$refGenotype)
    controlSubset <- GenotypeSubset[,c(depVariable)]
    controlSubset <- na.omit(controlSubset)
    
    # Reference interval is defined by using percentiles
    if (length(controlSubset)>=50){
        percentiles<-unname(quantile(controlSubset, c(.025, .975)))
        RR_left_all<-percentiles[1]
        RR_right_all<-percentiles[2]
    }
    else {
    # Reference interval is defined by suing Gaussian distribution
        # shapiro.test(controlSubset) is needed?    
        RR_left_all <- mean(controlSubset) - sd(controlSubset)
        RR_right_all <- mean(controlSubset) + sd(controlSubset)
    }
    
    if(numberofgenders==2){        
        
        GenotypeSubset_male <- subset(GenotypeSubset, 
                GenotypeSubset$Gender=="Male")
        
        GenotypeSubset_female <- subset(GenotypeSubset, 
                GenotypeSubset$Gender=="Female")
        
        controlSubset_male <- GenotypeSubset_male[,c(depVariable)]
        controlSubset_male <- na.omit(controlSubset_male)
        
        # Reference interval is defined by using percentiles
        if (length(controlSubset_male)>=50){
            percentiles<-unname(quantile(controlSubset_male, c(.025, .975)))
            RR_left_male<-percentiles[1]
            RR_right_male<-percentiles[2]
        }
        else {
            # Reference interval is defined by suing Gaussian distribution
            # shapiro.test(controlSubset) is needed?    
            RR_left_male <- mean(controlSubset_male) - sd(controlSubset_male)
            RR_right_male <- mean(controlSubset_male) + sd(controlSubset_male)
        }
        
        controlSubset_female <- GenotypeSubset_female[,c(depVariable)]
        controlSubset_female <- na.omit(controlSubset_female)
        
        # Reference interval is defined by using percentiles
        if (length(controlSubset_female)>=50){
            percentiles<-unname(quantile(controlSubset_female, c(.025, .975)))
            RR_left_female<-percentiles[1]
            RR_right_female<-percentiles[2]
        }
        else {
            # Reference interval is defined by suing Gaussian distribution
            # shapiro.test(controlSubset) is needed?    
            RR_left_female <- mean(controlSubset_female) - sd(controlSubset_female)
            RR_right_female <- mean(controlSubset_female) + sd(controlSubset_female)
        }
        
    }
    
    # RR plus test - all together
    GenotypeSubset <- subset(x, x$Genotype==phenList$testGenotype)
    mutantSubset <- GenotypeSubset[,c(depVariable)]
    mutantSubset <- na.omit(mutantSubset)
    valuesInRange <- sum(sapply(mutantSubset, function(x) (x>=RR_left && x<=RR_right)))
    allValues <- length(mutantSubset)
    percentageOut_all<-round(100*(allValues-valuesInRange)/allValues)
    
 
    
    
    # RR plus test - males/females
    if(numberofgenders==2){       
        
        GenotypeSubset_male <- subset(GenotypeSubset, 
                GenotypeSubset$Gender=="Male")
        
        GenotypeSubset_female <- subset(GenotypeSubset, 
                GenotypeSubset$Gender=="Female")
        
        mutantSubset_male <- GenotypeSubset_male[,c(depVariable)]
        mutantSubset_male <- na.omit(mutantSubset_male)
        
        mutantSubset_female <- GenotypeSubset_female[,c(depVariable)]
        mutantSubset_female <- na.omit(mutantSubset_female)
        
        valuesInRange_male <- sum(sapply(mutantSubset_male, function(x) (x>=RR_left_male && x<=RR_right_male)))
        allValues_male <- length(mutantSubset_male)
        percentageOut_male<-round(100*(allValues_male-valuesInRange_male)/allValues_male)
        
        valuesInRange_female <- sum(sapply(mutantSubset_female, function(x) (x>=RR_left_female && x<=RR_right_female)))
        allValues_female <- length(mutantSubset_female)
        percentageOut_female<-round(100*(allValues_female-valuesInRange_female)/allValues_female)
    }
    
    
    
    stat_all_list <- c(RR_left_all,RR_right_all,valuesInRange,allValues,percentageOut_all)
    #names(stat_all_list) <- c("RR left", "RR right"," values in range", "all values", "out of range %")
    
    stat_male_list <- c(RR_left_male,RR_right_male,valuesInRange_male,allValues_male,percentageOut_male)
    #names(stat_male_list) <- c("RR left", "RR right"," values in range", "all values", "out of range %")
    
    stat_female_list <- c(RR_left_female,RR_right_female,valuesInRange_female,allValues_female,percentageOut_female)
    #names(stat_female_list) <- c("RR left", "RR right"," values in range", "all values", "out of range %")
    
    
    estimateAll <- "Genotype effect is significant"
    if (percentageOut_all<rangeThreshold)
        estimateAll <- "Genotype effect is not significant"
    
    estimateMale <- "Genotype effect is significant for males only"
    if (percentageOut_male<rangeThreshold)
        estimateMale <- "Genotype effect is not significant for males only"
    
    estimateFemale <- "Genotype effect is significant for females only"
    if (percentageOut_female<rangeThreshold)
        estimateFemale <- "Genotype effect is not significant for females only"
    

    model <- c(estimateAll,  stat_all_list, estimateMale, stat_male_list, estimateFemale, stat_female_list)
    
    names(model) <- c("RR estimate", "RR left", "RR right", "values in range", "all values", "out of range %", 
            "RR estimate males only", "RR left males only", "RR right males only"," values in range males only", 
            "all values males only", "out of range % males only",
            "RR estimate females only", "RR left females only", "RR right females only"," values in range females only", 
            "all values females only", "out of range % females only")

    
    keep_weight <- NA
    keep_gender <- NA
    keep_interaction <- NA
    keep_batch <- NA
    keep_equalvar <- NA
    interactionTest <- NA
    

    result <- new("PhenTestResult",list(model.dataset=x, model.output=model,
                    depVariable=depVariable,method="RR",model.effect.batch=keep_batch,
                    model.effect.variance=keep_equalvar,model.effect.interaction=keep_interaction,
                    model.output.interaction=interactionTest,model.effect.gender=keep_gender,
                    model.effect.weight=keep_weight,numberGenders=numberofgenders))
    return(result)
}



##------------------------------------------------------------------------------
## Parser model output summary and return in readable vector format
parserOutputSummaryRR<-function(phenTestResult)
{
    result <- phenTestResult
    model <- result$model.output
    estimate <-model[1]
    RR_left <- model[2]
    RR_right <- model[3]
    values_in_range <- model[4]
    all_values <- model[5]
    values_out_percentage <- model[6]    
    
   
    output <- c(estimate,RR_left, RR_right, values_in_range, all_values, values_out_percentage)
    
    names(output) <- c("RR estimate","RR left", "RR right", "values in range", "all values", "out of range %")
    
    return(output)
}
##------------------------------------------------------------------------------