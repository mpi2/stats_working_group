## Copyright © 2013-2014 EMBL - European Bioinformatics Institute
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
# Non-Gaussian distribution - percentiles are used to limit of the reference range. 
# Sufficient numbers of control samples are required (minimum 60) - tests in testDataset function
# Natural variation to default to 95% min 75% and max 99%

RRTest <- function(phenList, depVariable, 
        outputMessages=TRUE, naturalVariation=95, controlPointsThreshold=60)
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
    rangesVector <- NA
    
    numberofsexes <- length(levels(x$Sex))
    
    # Define control subset
    GenotypeControlSubset <- subset(x, x$Genotype==phenList$refGenotype)
    
    
    # Define mutant subset
    GenotypeMutantSubset <- subset(x, x$Genotype==phenList$testGenotype)
    
    # Calculate reference ranges
    rangeLeft <- (100-naturalVariation)/200 
    rangeRight <- (100-((100-naturalVariation)/2))/100
    
    #allValues_male <- length(mutantSubset_male)
    if(numberofsexes==2){        
        # Reference interval for males
        GenotypeSubset_males <- subset(GenotypeControlSubset, 
                GenotypeControlSubset$Sex=="Male")
        controlSubset_males <- GenotypeSubset_males[,c(depVariable)]
        controlSubset_males <- na.omit(controlSubset_males)
        percentiles<-unname(quantile(controlSubset_males, c(rangeLeft, rangeRight)))
        RR_left_male<-percentiles[1]
        RR_right_male<-percentiles[2]
        
        countsNormal_control_males <- sum(sapply(controlSubset_males, 
                        function(x) (x>RR_left_male && x<RR_right_male)))
        countsLow_control_males <- sum(sapply(controlSubset_males, 
                        function(x) (x<=RR_left_male)))
        countsHigh_control_males <- sum(sapply(controlSubset_males, 
                        function(x) (x>=RR_right_male)))
        
        # Counts for males
        GenotypeSubset_males <- subset(GenotypeMutantSubset, 
                GenotypeMutantSubset$Sex=="Male")                      
        mutantSubset_males <- GenotypeSubset_males[,c(depVariable)]
        mutantSubset_males <- na.omit(mutantSubset_males)            
        countsNormal_mutant_males <- sum(sapply(mutantSubset_males, 
                        function(x) (x>RR_left_male && x<RR_right_male)))
        countsLow_mutant_males <- sum(sapply(mutantSubset_males, 
                        function(x) (x<=RR_left_male)))
        countsHigh_mutant_males <- sum(sapply(mutantSubset_males, 
                        function(x) (x>=RR_right_male)))
        
        # Reference interval for females
        GenotypeSubset_females <- subset(GenotypeControlSubset, 
                GenotypeControlSubset$Sex=="Female")
        controlSubset_females <- GenotypeSubset_females[,c(depVariable)]
        controlSubset_females <- na.omit(controlSubset_females)
        percentiles<-unname(quantile(controlSubset_females, c(rangeLeft, rangeRight)))
        RR_left_female<-percentiles[1]
        RR_right_female<-percentiles[2]
        
        countsNormal_control_females <- sum(sapply(controlSubset_females, 
                        function(x) (x>RR_left_female && x<RR_right_female)))
        countsLow_control_females <- sum(sapply(controlSubset_females, 
                        function(x) (x<=RR_left_female)))
        countsHigh_control_females <- sum(sapply(controlSubset_females, 
                        function(x) (x>=RR_right_female)))
        
        
        # Counts for females
        GenotypeSubset_females <- subset(GenotypeMutantSubset, 
                GenotypeMutantSubset$Sex=="Female")            
        mutantSubset_females <- GenotypeSubset_females[,c(depVariable)]
        mutantSubset_females <- na.omit(mutantSubset_females)            
        countsNormal_mutant_females <- sum(sapply(mutantSubset_females, 
                        function(x) (x>RR_left_female && x<RR_right_female)))
        countsLow_mutant_females <- sum(sapply(mutantSubset_females, 
                        function(x) (x<=RR_left_female)))
        countsHigh_mutant_females <- sum(sapply(mutantSubset_females, 
                        function(x) (x>=RR_right_female)))
        
        # Counts for all are sums of males and females values
        countsNormal_control_all <- countsNormal_control_males + countsNormal_control_females
        countsLow_control_all <- countsLow_control_males + countsLow_control_females
        countsHigh_control_all <- countsHigh_control_males + countsHigh_control_females 
        countsNormal_mutant_all <- countsNormal_mutant_males + countsNormal_mutant_females
        countsLow_mutant_all <- countsLow_mutant_males + countsLow_mutant_females
        countsHigh_mutant_all <- countsHigh_mutant_males + countsHigh_mutant_females
        
        rangesVector <- c(RR_left_male,RR_right_male,RR_left_female,RR_right_female)
    }
    
    else {
        # Reference interval for all
        controlSubset <- GenotypeControlSubset[,c(depVariable)]
        controlSubset <- na.omit(controlSubset)
        percentiles<-unname(quantile(controlSubset, c(rangeLeft, rangeRight)))
        RR_left_all<-percentiles[1]
        RR_right_all<-percentiles[2]
        rangesVector <- c(RR_left_all,RR_right_all)
        
        countsNormal_control_all <- sum(sapply(controlSubset, function(x) (x>RR_left_all && x<RR_right_all)))
        countsLow_control_all <- sum(sapply(controlSubset, function(x) (x<=RR_left_all)))
        countsHigh_control_all <- sum(sapply(controlSubset, function(x) (x>=RR_right_all)))
        
        
        # Counts for all
        mutantSubset <- GenotypeMutantSubset[,c(depVariable)]
        mutantSubset <- na.omit(mutantSubset)
        countsNormal_mutant_all <- sum(sapply(mutantSubset, function(x) (x>RR_left_all && x<RR_right_all)))
        countsLow_mutant_all <- sum(sapply(mutantSubset, function(x) (x<=RR_left_all)))
        countsHigh_mutant_all <- sum(sapply(mutantSubset, function(x) (x>=RR_right_all)))
        
    }
    
    # Matrices for combined datasets or one sex only
    
    #ALL
    
    #1) ##### ALL three rows: Low, Normal, High #####
    count_matrix_all <- matrix(0,3,2)
    ES_matrix_all <- matrix(0,3,3)
    count_matrix_all[1,1] <- countsLow_control_all
    count_matrix_all[2,1] <- countsNormal_control_all
    count_matrix_all[3,1] <- countsHigh_control_all
    count_matrix_all[1,2] <- countsLow_mutant_all
    count_matrix_all[2,2] <- countsNormal_mutant_all
    count_matrix_all[3,2] <- countsHigh_mutant_all
    ES_matrix_all[1,1] <- round((count_matrix_all[1,1]/colSums(count_matrix_all)[1])*100,digits=0)
    ES_matrix_all[2,1] <- round((count_matrix_all[2,1]/colSums(count_matrix_all)[1])*100,digits=0)
    ES_matrix_all[3,1] <- round((count_matrix_all[3,1]/colSums(count_matrix_all)[1])*100,digits=0)    
    ES_matrix_all[1,2] <- round((count_matrix_all[1,2]/colSums(count_matrix_all)[2])*100,digits=0)
    ES_matrix_all[2,2] <- round((count_matrix_all[2,2]/colSums(count_matrix_all)[2])*100,digits=0)
    ES_matrix_all[3,2] <- round((count_matrix_all[3,2]/colSums(count_matrix_all)[2])*100,digits=0)
    ES_matrix_all[1,3] <- abs(ES_matrix_all[1,1]-ES_matrix_all[1,2])
    ES_matrix_all[2,3] <- abs(ES_matrix_all[2,1]-ES_matrix_all[2,2])
    ES_matrix_all[3,3] <- abs(ES_matrix_all[3,1]-ES_matrix_all[3,2])
    colnames(count_matrix_all) <- c(phenList$refGenotype,phenList$testGenotype)
    rownames(count_matrix_all) <- c("Low","Normal","High")
    colnames(ES_matrix_all) <- c(phenList$refGenotype,phenList$testGenotype,"ES change")
    rownames(ES_matrix_all) <- c("Low","Normal","High")
    ###### ALL three rows: Low, Normal, High ######
    
    #2) ##### ALL two rows: Low, Normal/High #####
    count_matrix_all_nh <- matrix(0,2,2)
    ES_matrix_all_nh <- matrix(0,2,3)
    count_matrix_all_nh[1,1] <- countsLow_control_all
    count_matrix_all_nh[2,1] <- countsNormal_control_all + countsHigh_control_all
    count_matrix_all_nh[1,2] <- countsLow_mutant_all
    count_matrix_all_nh[2,2] <- countsNormal_mutant_all + countsHigh_mutant_all
    ES_matrix_all_nh[1,1] <- round((count_matrix_all_nh[1,1]/colSums(count_matrix_all_nh)[1])*100,digits=0)
    ES_matrix_all_nh[2,1] <- round((count_matrix_all_nh[2,1]/colSums(count_matrix_all_nh)[1])*100,digits=0)
    ES_matrix_all_nh[1,2] <- round((count_matrix_all_nh[1,2]/colSums(count_matrix_all_nh)[2])*100,digits=0)
    ES_matrix_all_nh[2,2] <- round((count_matrix_all_nh[2,2]/colSums(count_matrix_all_nh)[2])*100,digits=0)
    ES_matrix_all_nh[1,3] <- abs(ES_matrix_all_nh[1,1]-ES_matrix_all_nh[1,2])
    ES_matrix_all_nh[2,3] <- abs(ES_matrix_all_nh[2,1]-ES_matrix_all_nh[2,2])
    colnames(count_matrix_all_nh) <- c(phenList$refGenotype,phenList$testGenotype)
    rownames(count_matrix_all_nh) <- c("Low","Normal/High")
    colnames(ES_matrix_all_nh) <- c(phenList$refGenotype,phenList$testGenotype,"ES change")
    rownames(ES_matrix_all_nh) <- c("Low","Normal/High")
    ES_all_nh <- round(ES_matrix_all_nh[1,3],digits=0)
    model_all_nh <- fisher.test(count_matrix_all_nh)
    ###### ALL two rows: Low, Normal/High ######
    
    #3) ##### ALL two rows: High, Normal/Low #####
    count_matrix_all_nl <- matrix(0,2,2)
    ES_matrix_all_nl <- matrix(0,2,3)    
    count_matrix_all_nl[1,1] <- countsHigh_control_all
    count_matrix_all_nl[2,1] <- countsNormal_control_all + countsLow_control_all
    count_matrix_all_nl[1,2] <- countsHigh_mutant_all
    count_matrix_all_nl[2,2] <- countsNormal_mutant_all + countsLow_mutant_all
    ES_matrix_all_nl[1,1] <- round((count_matrix_all_nl[1,1]/colSums(count_matrix_all_nl)[1])*100,digits=0)
    ES_matrix_all_nl[2,1] <- round((count_matrix_all_nl[2,1]/colSums(count_matrix_all_nl)[1])*100,digits=0)
    ES_matrix_all_nl[1,2] <- round((count_matrix_all_nl[1,2]/colSums(count_matrix_all_nl)[2])*100,digits=0)
    ES_matrix_all_nl[2,2] <- round((count_matrix_all_nl[2,2]/colSums(count_matrix_all_nl)[2])*100,digits=0)
    ES_matrix_all_nl[1,3] <- abs(ES_matrix_all_nl[1,1]-ES_matrix_all_nl[1,2])
    ES_matrix_all_nl[2,3] <- abs(ES_matrix_all_nl[2,1]-ES_matrix_all_nl[2,2])
    colnames(count_matrix_all_nl) <- c(phenList$refGenotype,phenList$testGenotype)
    rownames(count_matrix_all_nl) <- c("High","Normal/Low")
    colnames(ES_matrix_all_nl) <- c(phenList$refGenotype,phenList$testGenotype,"ES change")
    rownames(ES_matrix_all_nl) <- c("High","Normal/Low")
    # Statistics
    ES_all_nl <- round(ES_matrix_all_nl[1,3],digits=0)
    model_all_nl <- fisher.test(count_matrix_all_nl)
    ###### ALL two rows: High, Normal/Low ######    
    
    model <- NULL
    model_male <- NULL
    model_female <- NULL
    count_matrix_female <-NULL
    count_matrix_male <- NULL
    ES_male <- NULL
    ES_female <- NULL
    ES_matrix_male<-NULL
    ES_matrix_female<-NULL
    ES_male_nh <- NULL
    ES_male_nl <- NULL
    ES_female_nh <- NULL
    ES_female_nl <- NULL
    stat_male <-NULL
    stat_female <- NULL
    
    # MALE/FEMALE ONLY
    if(numberofsexes==2){     
        #1) ##### FEMALE/MALE ONLY three rows: Low, Normal, High #####  
        count_matrix_male <- matrix(0,3,2)
        count_matrix_female <- matrix(0,3,2)
        ES_matrix_male <- matrix(0,3,3)
        ES_matrix_female <- matrix(0,3,3) 
        count_matrix_male[1,1] <- countsLow_control_males
        count_matrix_male[2,1] <- countsNormal_control_males
        count_matrix_male[3,1] <- countsHigh_control_males
        count_matrix_male[1,2] <- countsLow_mutant_males
        count_matrix_male[2,2] <- countsNormal_mutant_males
        count_matrix_male[3,2] <- countsHigh_mutant_males
        count_matrix_female[1,1] <- countsLow_control_females
        count_matrix_female[2,1] <- countsNormal_control_females
        count_matrix_female[3,1] <- countsHigh_control_females
        count_matrix_female[1,2] <- countsLow_mutant_females
        count_matrix_female[2,2] <- countsNormal_mutant_females
        count_matrix_female[3,2] <- countsHigh_mutant_females
        ES_matrix_male[1,1] <- round((count_matrix_male[1,1]/colSums(count_matrix_male)[1])*100,digits=0)
        ES_matrix_male[2,1] <- round((count_matrix_male[2,1]/colSums(count_matrix_male)[1])*100,digits=0)
        ES_matrix_male[3,1] <- round((count_matrix_male[3,1]/colSums(count_matrix_male)[1])*100,digits=0)    
        ES_matrix_male[1,2] <- round((count_matrix_male[1,2]/colSums(count_matrix_male)[2])*100,digits=0)
        ES_matrix_male[2,2] <- round((count_matrix_male[2,2]/colSums(count_matrix_male)[2])*100,digits=0)
        ES_matrix_male[3,2] <- round((count_matrix_male[3,2]/colSums(count_matrix_male)[2])*100,digits=0)
        ES_matrix_male[1,3] <- abs(ES_matrix_male[1,1]-ES_matrix_male[1,2])
        ES_matrix_male[2,3] <- abs(ES_matrix_male[2,1]-ES_matrix_male[2,2])
        ES_matrix_male[3,3] <- abs(ES_matrix_male[3,1]-ES_matrix_male[3,2])
        ES_matrix_female[1,1] <- round((count_matrix_female[1,1]/colSums(count_matrix_female)[1])*100,digits=0)
        ES_matrix_female[2,1] <- round((count_matrix_female[2,1]/colSums(count_matrix_female)[1])*100,digits=0)
        ES_matrix_female[3,1] <- round((count_matrix_female[3,1]/colSums(count_matrix_female)[1])*100,digits=0)    
        ES_matrix_female[1,2] <- round((count_matrix_female[1,2]/colSums(count_matrix_female)[2])*100,digits=0)
        ES_matrix_female[2,2] <- round((count_matrix_female[2,2]/colSums(count_matrix_female)[2])*100,digits=0)
        ES_matrix_female[3,2] <- round((count_matrix_female[3,2]/colSums(count_matrix_female)[2])*100,digits=0)
        ES_matrix_female[1,3] <- abs(ES_matrix_female[1,1]-ES_matrix_female[1,2])
        ES_matrix_female[2,3] <- abs(ES_matrix_female[2,1]-ES_matrix_female[2,2])
        ES_matrix_female[3,3] <- abs(ES_matrix_female[3,1]-ES_matrix_female[3,2])
        colnames(count_matrix_male) <- c(phenList$refGenotype,phenList$testGenotype)
        rownames(count_matrix_male) <- c("Low","Normal","High")
        colnames(count_matrix_female) <- c(phenList$refGenotype,phenList$testGenotype)
        rownames(count_matrix_female) <- c("Low","Normal","High")
        colnames(ES_matrix_male) <- c(phenList$refGenotype,phenList$testGenotype,"ES change")
        rownames(ES_matrix_male) <- c("Low","Normal","High")
        colnames(ES_matrix_female) <- c(phenList$refGenotype,phenList$testGenotype,"ES change")
        rownames(ES_matrix_female) <- c("Low","Normal","High")
        ###### FEMALE/MALE ONLY three rows: Low, Normal, High ######
        
        #2) ##### FEMALE/MALE ONLY two rows: Low, Normal/High #####
        count_matrix_female_nh <- matrix(0,2,2)
        ES_matrix_female_nh <- matrix(0,2,3)
        count_matrix_male_nh <- matrix(0,2,2)
        ES_matrix_male_nh <- matrix(0,2,3)
        count_matrix_male_nh[1,1] <- countsLow_control_males
        count_matrix_male_nh[2,1] <- countsNormal_control_males + countsHigh_control_males
        count_matrix_male_nh[1,2] <- countsLow_mutant_males
        count_matrix_male_nh[2,2] <- countsNormal_mutant_males + countsHigh_mutant_males
        count_matrix_female_nh[1,1] <- countsLow_control_females
        count_matrix_female_nh[2,1] <- countsNormal_control_females + countsHigh_control_females
        count_matrix_female_nh[1,2] <- countsLow_mutant_females
        count_matrix_female_nh[2,2] <- countsNormal_mutant_females + countsHigh_mutant_females
        ES_matrix_male_nh[1,1] <- round((count_matrix_male_nh[1,1]/colSums(count_matrix_male_nh)[1])*100,digits=0)
        ES_matrix_male_nh[2,1] <- round((count_matrix_male_nh[2,1]/colSums(count_matrix_male_nh)[1])*100,digits=0)
        ES_matrix_male_nh[1,2] <- round((count_matrix_male_nh[1,2]/colSums(count_matrix_male_nh)[2])*100,digits=0)
        ES_matrix_male_nh[2,2] <- round((count_matrix_male_nh[2,2]/colSums(count_matrix_male_nh)[2])*100,digits=0)
        ES_matrix_male_nh[1,3] <- abs(ES_matrix_male_nh[1,1]-ES_matrix_male_nh[1,2])
        ES_matrix_male_nh[2,3] <- abs(ES_matrix_male_nh[2,1]-ES_matrix_male_nh[2,2])
        ES_matrix_female_nh[1,1] <- round((count_matrix_female_nh[1,1]/colSums(count_matrix_female_nh)[1])*100,digits=0)
        ES_matrix_female_nh[2,1] <- round((count_matrix_female_nh[2,1]/colSums(count_matrix_female_nh)[1])*100,digits=0)
        ES_matrix_female_nh[1,2] <- round((count_matrix_female_nh[1,2]/colSums(count_matrix_female_nh)[2])*100,digits=0)
        ES_matrix_female_nh[2,2] <- round((count_matrix_female_nh[2,2]/colSums(count_matrix_female_nh)[2])*100,digits=0)
        ES_matrix_female_nh[1,3] <- abs(ES_matrix_female_nh[1,1]-ES_matrix_female_nh[1,2])
        ES_matrix_female_nh[2,3] <- abs(ES_matrix_female_nh[2,1]-ES_matrix_female_nh[2,2])
        colnames(count_matrix_male_nh) <- c(phenList$refGenotype,phenList$testGenotype)
        rownames(count_matrix_male_nh) <- c("Low","Normal/High")
        colnames(ES_matrix_male_nh) <- c(phenList$refGenotype,phenList$testGenotype,"ES change")
        rownames(ES_matrix_male_nh) <- c("Low","Normal/High")
        colnames(count_matrix_female_nh) <- c(phenList$refGenotype,phenList$testGenotype)
        rownames(count_matrix_female_nh) <- c("Low","Normal/High")
        colnames(ES_matrix_female_nh) <- c(phenList$refGenotype,phenList$testGenotype,"ES change")
        rownames(ES_matrix_female_nh) <- c("Low","Normal/High")
        # Statistics
        ES_male_nh <- round(ES_matrix_male_nh[1,3],digits=0)
        model_male_nh <- fisher.test(count_matrix_male_nh)
        ES_female_nh <- round(ES_matrix_female_nh[1,3],digits=0)
        model_female_nh <- fisher.test(count_matrix_female_nh)
        ###### FEMALE/MALE ONLY two rows: Low, Normal/High ######
        
        #3) ##### FEMALE/MALE ONLY two rows: High, Normal/Low #####
        count_matrix_female_nl <- matrix(0,2,2)
        ES_matrix_female_nl <- matrix(0,2,3)
        count_matrix_male_nl <- matrix(0,2,2)
        ES_matrix_male_nl <- matrix(0,2,3)
        count_matrix_male_nl[1,1] <- countsHigh_control_males
        count_matrix_male_nl[2,1] <- countsNormal_control_males + countsLow_control_males
        count_matrix_male_nl[1,2] <- countsHigh_mutant_males
        count_matrix_male_nl[2,2] <- countsNormal_mutant_males + countsLow_mutant_males
        count_matrix_female_nl[1,1] <- countsHigh_control_females
        count_matrix_female_nl[2,1] <- countsNormal_control_females + countsLow_control_females
        count_matrix_female_nl[1,2] <- countsHigh_mutant_females
        count_matrix_female_nl[2,2] <- countsNormal_mutant_females + countsLow_mutant_females
        ES_matrix_male_nl[1,1] <- round((count_matrix_male_nl[1,1]/colSums(count_matrix_male_nl)[1])*100,digits=0)
        ES_matrix_male_nl[2,1] <- round((count_matrix_male_nl[2,1]/colSums(count_matrix_male_nl)[1])*100,digits=0)
        ES_matrix_male_nl[1,2] <- round((count_matrix_male_nl[1,2]/colSums(count_matrix_male_nl)[2])*100,digits=0)
        ES_matrix_male_nl[2,2] <- round((count_matrix_male_nl[2,2]/colSums(count_matrix_male_nl)[2])*100,digits=0)
        ES_matrix_male_nl[1,3] <- abs(ES_matrix_male_nl[1,1]-ES_matrix_male_nl[1,2])
        ES_matrix_male_nl[2,3] <- abs(ES_matrix_male_nl[2,1]-ES_matrix_male_nl[2,2])
        ES_matrix_female_nl[1,1] <- round((count_matrix_female_nl[1,1]/colSums(count_matrix_female_nl)[1])*100,digits=0)
        ES_matrix_female_nl[2,1] <- round((count_matrix_female_nl[2,1]/colSums(count_matrix_female_nl)[1])*100,digits=0)
        ES_matrix_female_nl[1,2] <- round((count_matrix_female_nl[1,2]/colSums(count_matrix_female_nl)[2])*100,digits=0)
        ES_matrix_female_nl[2,2] <- round((count_matrix_female_nl[2,2]/colSums(count_matrix_female_nl)[2])*100,digits=0)
        ES_matrix_female_nl[1,3] <- abs(ES_matrix_female_nl[1,1]-ES_matrix_female_nl[1,2])
        ES_matrix_female_nl[2,3] <- abs(ES_matrix_female_nl[2,1]-ES_matrix_female_nl[2,2])
        colnames(count_matrix_male_nl) <- c(phenList$refGenotype,phenList$testGenotype)
        rownames(count_matrix_male_nl) <- c("High","Normal/Low")
        colnames(ES_matrix_male_nl) <- c(phenList$refGenotype,phenList$testGenotype,"ES change")
        rownames(ES_matrix_male_nl) <- c("High","Normal/Low")
        colnames(count_matrix_female_nl) <- c(phenList$refGenotype,phenList$testGenotype)
        rownames(count_matrix_female_nl) <- c("High","Normal/Low")
        colnames(ES_matrix_female_nl) <- c(phenList$refGenotype,phenList$testGenotype,"ES change")
        rownames(ES_matrix_female_nl) <- c("High","Normal/Low")
        # Statistics
        ES_male_nl <- round(ES_matrix_male_nl[1,3],digits=0)
        model_male_nl <- fisher.test(count_matrix_male_nl)
        ES_female_nl <- round(ES_matrix_female_nl[1,3],digits=0)
        model_female_nl <- fisher.test(count_matrix_female_nl)
        ###### FEMALE/MALE ONLY two rows: High, Normal/Low ######
        
        #stat_male <- assocstats(count_matrix_male)
        #stat_female <- assocstats(count_matrix_female)
    }
    
    
    # Combine the results
    results_all <- matrix(0,4,1)
    results_all[1] <- sprintf("%.4f",round(p.adjust(model_all_nh$p.value, method = "bonferroni", n = 2),4))
    results_all[2] <- paste(format(ES_all_nh, nsmall = 0),"%",sep="")
    results_all[3] <- sprintf("%.4f",round(p.adjust(model_all_nl$p.value, method = "bonferroni", n = 2),4))
    results_all[4] <- paste(format(ES_all_nl, nsmall = 0),"%",sep="")
    rownames(results_all) <- c("Low classification p-value:","Low classification effect size:",
            "High classification p-value:","High classification effect size:")
    colnames(results_all) <- c("")
    
    model$count_matrix_female <- count_matrix_female
    model$count_matrix_male <- count_matrix_male
    model$count_matrix_all <- count_matrix_all
    #model$stat_all <- assocstats(count_matrix_all)
    model$ES <- c(ES_all_nh,ES_all_nl)
    model$ES_male <- c(ES_male_nh, ES_male_nl)
    model$ES_female <- c(ES_female_nh,ES_female_nl)
    model$percentage_matrix_all <-ES_matrix_all
    model$percentage_matrix_male <-ES_matrix_male
    model$percentage_matrix_female <-ES_matrix_female
    #model$stat_male <- stat_male
    #model$stat_female <- stat_female
    
    if(numberofsexes==2){    
        thresholds <- matrix(0,4,1)
        thresholds[1] <- format(naturalVariation, nsmall = 0)
        thresholds[2] <- format(controlPointsThreshold, nsmall = 0)
        thresholds[3] <- paste(format(RR_left_male,nsmall=3)," to ",format(RR_right_male,nsmall=3),sep="")
        thresholds[4] <- paste(format(RR_left_female,nsmall=3)," to ",format(RR_right_female,nsmall=3),sep="")
        rownames(thresholds) <- c("Natural variation:","Min control points:",
                "Normal values 'males only':","Normal values 'females only':")
        colnames(thresholds) <- c("")
        
        results_male <- matrix(0,4,1)
        results_male[1] <- sprintf("%.4f",round(p.adjust(model_male_nh$p.value, method = "bonferroni", n = 2),4))
        results_male[2] <- paste(format(ES_male_nh, nsmall = 0),"%",sep="")
        results_male[3] <- sprintf("%.4f",round(p.adjust(model_male_nl$p.value, method = "bonferroni", n = 2),4))
        results_male[4] <- paste(format(ES_male_nl, nsmall = 0),"%",sep="")
        rownames(results_male) <- c("Low classification p-value:","Low classification effect size:",
                "High classification p-value:","High classification effect size:")
        colnames(results_male) <- c("")
        
        results_female <- matrix(0,4,1)
        results_female[1] <- sprintf("%.4f",round(p.adjust(model_female_nh$p.value, method = "bonferroni", n = 2),4))
        results_female[2] <- paste(format(ES_female_nh, nsmall = 0),"%",sep="")
        results_female[3] <- sprintf("%.4f",round(p.adjust(model_female_nl$p.value, method = "bonferroni", n = 2),4))
        results_female[4] <- paste(format(ES_female_nl, nsmall = 0),"%",sep="")
        rownames(results_female) <- c("Low classification p-value:","Low classification effect size:",
                "High classification p-value:","High classification effect size:")
        colnames(results_female) <- c("")
    }
    else {
        thresholds <- matrix(0,3,1)
        thresholds[1] <- format(naturalVariation, nsmall = 0)
        thresholds[2] <- format(controlPointsThreshold, nsmall = 0)
        thresholds[3] <- paste(format(RR_left_all,nsmall=3)," to ",format(RR_right_all,nsmall=3),sep="")
        rownames(thresholds) <- c("Natural variation:","Min control points:","Normal values 'all':")
        colnames(thresholds) <- c("")
        
        results_male <- NULL
        results_female <- NULL
    }
    
    model$all <- results_all
    model$male <- results_male
    model$female <- results_female
    
    keep_weight <- NA
    keep_sex <- NA
    keep_interaction <- NA
    keep_batch <- NA
    keep_equalvar <- NA
    interactionTest <- NA
    
    
    result <- new("PhenTestResult",list(model.dataset=x, model.output=model,
                    depVariable=depVariable,method="RR",model.effect.batch=keep_batch,
                    model.effect.variance=keep_equalvar,model.effect.interaction=keep_interaction,
                    model.output.interaction=interactionTest,model.effect.sex=keep_sex,
                    model.effect.weight=keep_weight,numberSexes=numberofsexes, model.output.quality = thresholds ))
    return(result)
}

##------------------------------------------------------------------------------