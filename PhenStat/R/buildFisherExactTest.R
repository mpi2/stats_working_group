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
# buildFisherExactTest.R contains buildFisherExactTest function
#-----------------------------------------------------------------------------------
# Create two times n table with record counts, where two rows represent genotype levels; n columns 
# represent dependent variable levels. Perform Fisher Exact Test
buildFisherExactTest <- function(phenList, depVariable, outputMessages=TRUE)
{    

    x <- phenList$dataset
    numberofgenders=length(levels(x$Gender))
    Genotype_levels=levels(x$Genotype)         
    depVariable_levels=levels(factor(x[,c(depVariable)]))  
    count_matrix_all <- matrix(0,2,length(depVariable_levels))
    
    for (i in 1:length(Genotype_levels)){
        GenotypeSubset <- subset(x, x$Genotype==Genotype_levels[i])
        for (j in 1:length(depVariable_levels)){  
            columnOfInterest <- GenotypeSubset[,c(depVariable)]
            columnOfInterest <- columnOfInterest[columnOfInterest==depVariable_levels[j]]
            nr <- length(columnOfInterest) 
            if (is.null(nr)) nr<-0 
            count_matrix_all[i,j]=nr
        }    
        
    }
    
    model_all<-fisher.test(count_matrix_all)
    
    rownames(count_matrix_all)<-Genotype_levels
    colnames(count_matrix_all)<-depVariable_levels

    
    if (numberofgenders==2){
        count_matrix_male <- matrix(0,2,length(depVariable_levels))
        count_matrix_female <- matrix(0,2,length(depVariable_levels))
        
        for (i in 1:length(Genotype_levels)){
            GenotypeSubset <- subset(x, x$Genotype==Genotype_levels[i])
            GenotypeSubset_male <- subset(GenotypeSubset, GenotypeSubset$Gender=="Male")
            GenotypeSubset_female <- subset(GenotypeSubset, GenotypeSubset$Gender=="Female")
            for (j in 1:length(depVariable_levels)){  
                columnOfInterest_male <- GenotypeSubset_male[,c(depVariable)]
                columnOfInterest_male <- columnOfInterest_male[columnOfInterest_male==depVariable_levels[j]]
                m_nr <- length(columnOfInterest_male) 
                if (is.null(m_nr) || is.na(m_nr)) m_nr<-0 
                count_matrix_male[i,j]=m_nr
                
                columnOfInterest_female <- GenotypeSubset_female[,c(depVariable)]
                columnOfInterest_female <- columnOfInterest_female[columnOfInterest_female==depVariable_levels[j]]
                f_nr <- length(columnOfInterest_female) 
                if (is.null(f_nr)) f_nr<-0 
                count_matrix_female[i,j]=f_nr
            }    
            
        } 
        
       
        rownames(count_matrix_female)<-Genotype_levels
        rownames(count_matrix_male)<-Genotype_levels
        colnames(count_matrix_female)<-depVariable_levels
        colnames(count_matrix_male)<-depVariable_levels
        model_male<-fisher.test(count_matrix_male)
        model_female<-fisher.test(count_matrix_female)
        model_male_results <- c(model_male$p.value,model_male$alternative,paste(model_male$conf.int[1],model_male$conf.int[2],sep=" to "),model_male$estimate)
        #names(model_male_results) <- c("Male p-value","Male alternative","Male confident interval","Male estimate")
        model_female_results <- c(model_female$p.value,model_female$alternative,paste(model_female$conf.int[1],model_female$conf.int[2],sep=" to "),model_female$estimate)  
        #names(model_female_results) <- c("Female p-value","Female alternative","Female confident interval","Female estimate")  
        
    }
    else {
        model_male_results <- c(NA,NA,NA,NA,NA)
        model_female_results <- c(NA,NA,NA,NA) 
        model_male <- NULL
        model_female <- NULL
        count_matrix_female <-NULL
        count_matrix_male <- NULL
    }    
    
    model_results <- c(model_all$p.value,model_all$alternative,paste(model_all$conf.int[1],model_all$conf.int[2],sep=" to "),model_all$estimate)
    #names(model_results) <- c("All p-value","All alternative","All confident interval","All estimate")
    
    #model$all <- model_results
    #model$male <- model_male_results
    #model$female <- model_female_results
    model <- model_all
    
    model$all <- model_all
    model$male <- model_male
    model$female <- model_female
    model$count_matrix_female <- count_matrix_female
    model$count_matrix_male <- count_matrix_male
    model$count_matrix_all <- count_matrix_all
    
    keep_weight <- NA
    keep_gender <- NA
    keep_interaction <- NA
    keep_batch <- NA
    keep_equalvar <- NA
    interactionTest <- NA
    
    result <- new("PhenTestResult",list(model.output=model,depVariable=depVariable,method="FE", 
                    model.effect.batch=keep_batch,model.effect.variance=keep_equalvar,model.effect.interaction=keep_interaction,
                    model.output.interaction=interactionTest,model.effect.gender=keep_gender,model.effect.weight=keep_weight,
                    numberGenders=numberofgenders))
    return(result)
}    