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
## FETFramework.R contains buildFisherExactTest function
##------------------------------------------------------------------------------
## Create two times n table with record counts, where two rows represent 
## genotype levels; n columns represent dependent variable levels. 
## Perform Fisher Exact test(s)
FisherExactTest <- function(phenList, depVariable, outputMessages=TRUE)
{    
    x <- getDataset(phenList)
    resultList <- list()
    indexList <- 1
    
    numberofsexes <- noSexes(phenList)
    Genotype_levels <- levels(x$Genotype)         
    depVariable_levels <- levels(factor(x[,c(depVariable)])) 
    
    if (length(depVariable_levels)==1){
        depVariable_levels <- c(depVariable_levels,"Other")
    }
    
    count_matrix_all <- matrix(0,length(depVariable_levels),2)
    ES_matrix_all <- matrix(0,length(depVariable_levels),3)
    
    
    for (i in 1:length(Genotype_levels)){
        GenotypeSubset <- subset(x, x$Genotype==Genotype_levels[i])
        for (j in 1:length(depVariable_levels)){  
            columnOfInterest <- GenotypeSubset[,c(depVariable)]
            columnOfInterest <- 
            columnOfInterest[columnOfInterest==depVariable_levels[j]]    
            
            columnOfInterest <- na.omit(columnOfInterest)
            nr <- length(columnOfInterest) 
            if (is.null(nr)) nr<-0 
            count_matrix_all[j,i] <- nr
        }    
        
    }
    
    for (i in 1:length(Genotype_levels)){
        GenotypeSubset <- subset(x, x$Genotype==Genotype_levels[i])
        for (j in 1:length(depVariable_levels)){  
            ES_matrix_all[j,i] <- round((count_matrix_all[j,i]/
                            colSums(count_matrix_all)[i])*100,digits=0)
        }    
        
    }
    for (j in 1:length(depVariable_levels)){  
        ES_matrix_all[j,3] <- abs(ES_matrix_all[j,1]-ES_matrix_all[j,2])
    }  
    
    colnames(count_matrix_all) <- Genotype_levels
    rownames(count_matrix_all) <- depVariable_levels
    
    colnames(ES_matrix_all) <- append(Genotype_levels,"ES change")
    rownames(ES_matrix_all) <- depVariable_levels
    
    ES_all <- round(max(ES_matrix_all[,3]),digits=0)
    
    model_all <- fisher.test(count_matrix_all)
    
    resultList[[indexList]] <- new("htestPhenStat",
            modelOutput=model_all,
            analysedSubset="all",
            comparison=character(0),
            matrixCount=count_matrix_all,
            ES=ES_all)
    indexList <- indexList+1
    
    model_male_results <- c(NA,NA,NA,NA,NA)
    model_female_results <- c(NA,NA,NA,NA) 
    model_male <- NULL
    model_female <- NULL
    count_matrix_female <-NULL
    count_matrix_male <- NULL
    ES_male <- NULL
    ES_female <- NULL
    ES_matrix_male<-NULL
    ES_matrix_female<-NULL
    stat_male <-NULL
    stat_female <- NULL
    
    
    if(numberofsexes==2){
        count_matrix_male <- matrix(0,length(depVariable_levels),2)
        count_matrix_female <- matrix(0,length(depVariable_levels),2)
        ES_matrix_male <- matrix(0,length(depVariable_levels),3)
        ES_matrix_female <- matrix(0,length(depVariable_levels),3)
        
        for (i in 1:length(Genotype_levels)){
            GenotypeSubset <- subset(x, x$Genotype==Genotype_levels[i])
            GenotypeSubset_male <- subset(GenotypeSubset, 
                    GenotypeSubset$Sex=="Male")
            
            GenotypeSubset_female <- subset(GenotypeSubset, 
                    GenotypeSubset$Sex=="Female")
            
            for (j in 1:length(depVariable_levels)){  
                columnOfInterest_male <- GenotypeSubset_male[,c(depVariable)]
                columnOfInterest_male <- 
                columnOfInterest_male[columnOfInterest_male==depVariable_levels[j]]
                
                columnOfInterest_male <- na.omit(columnOfInterest_male)
                m_nr <- length(columnOfInterest_male) 
                if (is.null(m_nr) || is.na(m_nr)) 
                m_nr <- 0 
                count_matrix_male[j,i] <- m_nr
                
                columnOfInterest_female <- GenotypeSubset_female[,c(depVariable)]
                columnOfInterest_female <- 
                columnOfInterest_female[columnOfInterest_female==depVariable_levels[j]]
                
                columnOfInterest_female <- na.omit(columnOfInterest_female)
                f_nr <- length(columnOfInterest_female) 
                if (is.null(f_nr)) 
                f_nr <- 0 
                count_matrix_female[j,i] <- f_nr
            }    
            
        } 
        
        
        for (i in 1:length(Genotype_levels)){
            GenotypeSubset <- subset(x, x$Genotype==Genotype_levels[i])
            for (j in 1:length(depVariable_levels)){  
                ES_matrix_male[j,i]=round((count_matrix_male[j,i]/
                                colSums(count_matrix_male)[i])*100,digits=0)
                ES_matrix_female[j,i]=round((count_matrix_female[j,i]/
                                colSums(count_matrix_female)[i])*100,digits=0)
            }    
            
        }
        for (j in 1:length(depVariable_levels)){  
            ES_matrix_male[j,3]=abs(ES_matrix_male[j,1]-ES_matrix_male[j,2])
            ES_matrix_female[j,3]=abs(ES_matrix_female[j,1]-ES_matrix_female[j,2])
        }  
        
        ES_male <- round(max(ES_matrix_male[,3]),digits=0)
        ES_female <- round(max(ES_matrix_female[,3]),digits=0)
        
        colnames(count_matrix_female) <- Genotype_levels
        colnames(count_matrix_male) <- Genotype_levels
        rownames(count_matrix_female) <- depVariable_levels
        rownames(count_matrix_male) <- depVariable_levels
        
        colnames(ES_matrix_male) <- append(Genotype_levels,"ES change")
        colnames(ES_matrix_female) <- append(Genotype_levels,"ES change")
        rownames(ES_matrix_male) <- depVariable_levels
        rownames(ES_matrix_female) <- depVariable_levels
        
        model_male <- fisher.test(count_matrix_male)
        model_female <- fisher.test(count_matrix_female)
        model_male_results <- c(model_male$p.value,model_male$alternative,
                paste(model_male$conf.int[1],model_male$conf.int[2],sep=" to "),
                model_male$estimate)
        
        model_female_results <- c(model_female$p.value,model_female$alternative,
                paste(model_female$conf.int[1],model_female$conf.int[2],sep=" to "),
                model_female$estimate)  
        
        #stat_male <- assocstats(count_matrix_male)
        #stat_female <- assocstats(count_matrix_female)
        
        resultList[[indexList]] <- new("htestPhenStat",
                modelOutput=model_male,
                analysedSubset="males",
                comparison=character(0),
                matrixCount=count_matrix_male,
                ES=ES_male)
        indexList <- indexList+1
        
        resultList[[indexList]] <- new("htestPhenStat",
                modelOutput=model_female,
                analysedSubset="females",
                comparison=character(0),
                matrixCount=count_matrix_female,
                ES=ES_female)
        indexList <- indexList+1
    }
    
    model_results <- c(model_all$p.value,model_all$alternative,
            paste(model_all$conf.int[1],model_all$conf.int[2],sep=" to "),
            model_all$estimate)
    
    model <- model_all
    
    model$all <- model_all
    model$male <- model_male
    model$female <- model_female
    model$count_matrix_female <- count_matrix_female
    model$count_matrix_male <- count_matrix_male
    model$count_matrix_all <- count_matrix_all
    #model$stat_all <- assocstats(count_matrix_all)
    model$ES <- ES_all
    model$ES_male <- ES_male
    model$ES_female <- ES_female
    model$percentage_matrix_all <-ES_matrix_all
    model$percentage_matrix_male <-ES_matrix_male
    model$percentage_matrix_female <-ES_matrix_female
    model$stat_male <- stat_male
    model$stat_female <- stat_female
    
    keep_weight <- NA
    keep_sex <- NA
    keep_interaction <- NA
    keep_batch <- NA
    keep_equalvar <- NA
    interactionTest <- NA
    
    
    #result <- new("PhenTestResult",list(model.dataset=x, model.output=model,
    #                depVariable=depVariable,refGenotype=phenList$refGenotype,method="FE",model.effect.batch=keep_batch,
    #                model.effect.variance=keep_equalvar,model.effect.interaction=keep_interaction,
    #                model.output.interaction=interactionTest,model.effect.sex=keep_sex,
    #                model.effect.weight=keep_weight,numberSexes=numberofsexes))
    

    result <- new("PhenTestResult", analysedDataset = x,
            depVariable=depVariable,
            refGenotype=refGenotype(phenList),
            testGenotype=testGenotype(phenList),
            method="FE",
            parameters=matrix(nrow=0, ncol=0),
            analysisResults = resultList)
    return(result)
} 
##------------------------------------------------------------------------------