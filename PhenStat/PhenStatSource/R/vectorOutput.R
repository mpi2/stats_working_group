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
## vectorOutput.R contains vectorOutput and outputLength functions
##------------------------------------------------------------------------------
## Wrapper to prepare the output of the modeling and testing results in vector 
## form. Assumes that modeling results are stored in the phenTestResult object 
## (output from functions testDataset and buildFinalModel)
vectorOutput <- function(phenTestResult)
{
 
    if (phenTestResult$method=="MM") {
        equation <- switch(phenTestResult$equation,
                withoutWeight = {"equation withoutWeight"},withWeight = {"equation withWeight"})
        
        fittingMethod <- "generalized least squares, "
        if (phenTestResult$model.effect.batch)
        fittingMethod <- "linear mixed-effects model, "
        
        classificationValue <- classificationTag(phenTestResult,
                userMode="vectorOutput",outputMessages=FALSE)
        
        
        x = phenTestResult$model.dataset
        columnOfInterest <- x[,c(phenTestResult$depVariable)]
        variability = paste('"variability":', round(length(unique(columnOfInterest))/length(columnOfInterest),digits=3),sep="")
        
        Genotype_levels=levels(x$Genotype)
        Sex_levels=levels(x$Sex)
        
        DSsize = ""
        
        for (i in 1:length(Genotype_levels)){
            GenotypeSubset <- subset(x, x$Genotype==Genotype_levels[i])
            for (j in 1:length(Sex_levels)){
                GenotypeSexSubset <- subset(GenotypeSubset, GenotypeSubset$Sex==Sex_levels[j])
                columnOfInterest <- GenotypeSexSubset[,c(phenTestResult$depVariable)]
                DSsize = paste(DSsize,paste('"',paste(Genotype_levels[i],Sex_levels[j],sep="_"),'"',sep=""),sep="")
                DSsize = paste(DSsize,":",sep="")
                DSsize = paste(DSsize,length(columnOfInterest),",",sep="")
            }
        }
        
        formula <- paste('"Formula":','"',
                        format(phenTestResult$model.formula.genotype),'"',sep="")
        
        addInfo = paste("{",DSsize,variability,",",formula,"}",sep="")
        
        vectorOutput <- c(paste("MM framework, ",fittingMethod, equation,sep=""),
                as.character(phenTestResult$depVariable), 
                as.character(phenTestResult$model.effect.batch), 
                as.character(phenTestResult$model.effect.variance), 
                as.character(phenTestResult$model.output.genotype.nulltest.pVal), 
                as.character(phenTestResult$model.output.summary["genotype_estimate"]), 
                as.character(phenTestResult$model.output.summary["genotype_estimate_SE"]),  
                as.character(phenTestResult$model.output.summary["genotype_p_value"]),
                as.character(phenTestResult$model.output.summary["sex_estimate"]), 
                as.character(phenTestResult$model.output.summary["sex_estimate_SE"]),  
                as.character(phenTestResult$model.output.summary["sex_p_value"]), 
                as.character(phenTestResult$model.output.summary["weight_estimate"]), 
                as.character(phenTestResult$model.output.summary["weight_estimate_SE"]), 
                as.character(phenTestResult$model.output.summary["weight_p_value"]), 
                as.character(phenTestResult$model.output.quality), 
                as.character(phenTestResult$model.output.summary["intercept_estimate"]), 
                as.character(phenTestResult$model.output.summary["intercept_estimate_SE"]), 
                as.character(phenTestResult$model.effect.interaction),
                as.character(phenTestResult$model.output.interaction),
                as.character(phenTestResult$model.output.summary["sex_FvKO_estimate"]), 
                as.character(phenTestResult$model.output.summary["sex_FvKO_SE"]), 
                as.character(phenTestResult$model.output.summary["sex_FvKO_p_value"]),  
                as.character(phenTestResult$model.output.summary["sex_MvKO_estimate"]),  
                as.character(phenTestResult$model.output.summary["sex_MvKO_SE"]), 
                as.character(phenTestResult$model.output.summary["sex_MvKO_p_value"]),
                as.character(classificationValue),
                as.character(addInfo))
        
        names(vectorOutput) <- c("Method",
                "Dependent variable",
                "Batch included",
                "Residual variances homogeneity",
                "Genotype contribution",
                "Genotype estimate",
                "Genotype standard error",
                "Genotype p-Val",
                "Sex estimate",
                "Sex standard error",
                "Sex p-val",
                "Weight estimate",
                "Weight standard error",
                "Weight p-val",
                "Gp1 genotype",
                "Gp1 Residuals normality test",
                "Gp2 genotype",
                "Gp2 Residuals normality test",
                "Blups test",
                "Rotated residuals normality test",
                "Intercept estimate",
                "Intercept standard error",
                "Interaction included",
                "Interaction p-val",
                "Sex FvKO estimate",
                "Sex FvKO standard error",
                "Sex FvKO p-val",
                "Sex MvKO estimate",
                "Sex MvKO standard error",
                "Sex MvKO p-val",
                "Classification tag",
                "Additional information")
    }
    else if (phenTestResult$method %in% c("FE")){
        male_pval <- NA
        female_pval <- NA
        male_ES <- NA
        female_ES <- NA
        if (!is.null(phenTestResult$model.output$male)){
            male_pval<-as.numeric(phenTestResult$model.output$male$p.val)
            male_ES<-as.numeric(phenTestResult$model.output$ES_male)
        }
        if (!is.null(phenTestResult$model.output$female)){
            female_pval<-as.numeric(phenTestResult$model.output$female$p.val)
            female_ES<-as.numeric(phenTestResult$model.output$ES_female)
        }
        
        classificationValue <- classificationTag(phenTestResult,
                userMode="vectorOutput",outputMessages=FALSE)
        
        vectorOutput <- c(switch(phenTestResult$method,FE = "Fisher Exact Test framework",
                        RR = "Reference Ranges Plus framework"),
                as.character(phenTestResult$depVariable), 
                "NA", 
                "NA",
                "NA",
                as.character(as.numeric(phenTestResult$model.output$ES)),
                "NA",  
                as.character(as.numeric(phenTestResult$model.output$all$p.val)),
                "NA",
                "NA", #10 
                "NA", 
                "NA", 
                "NA", 
                "NA", 
                as.character(colnames(phenTestResult$model.output$count_matrix_all)[1]), 
                "NA",#16 
                as.character(colnames(phenTestResult$model.output$count_matrix_all)[2]),  
                "NA", 
                "NA",
                "NA", 
                "NA", 
                "NA",
                "NA",
                "NA", #24
                as.character(female_ES),  
                "NA",  
                as.character(female_pval), 
                as.character(male_ES),  
                "NA",  
                as.character(male_pval), #30
                as.character(classificationValue),
                "NA")
        
        names(vectorOutput) <- c("Method",
                "Dependent variable",
                "Batch included",
                "Residual variances homogeneity",
                "Genotype contribution",
                "Genotype estimate",
                "Genotype standard error",
                "Genotype p-Val",
                "Sex estimate",
                "Sex standard error",
                "Sex p-val", #11
                "Weight estimate",
                "Weight standard error",
                "Weight p-val",
                "Gp1 genotype",
                "Gp1 Residuals normality test", #16
                "Gp2 genotype",
                "Gp2 Residuals normality test",
                "Blups test",
                "Rotated residuals normality test", #20
                "Intercept estimate",
                "Intercept standard error",
                "Interaction included",
                "Interaction p-val", #24
                "Sex FvKO estimate",
                "Sex FvKO standard error",
                "Sex FvKO p-val",
                "Sex MvKO estimate",
                "Sex MvKO standard error",
                "Sex MvKO p-val",
                "Classification tag",
                "Additional information")
    }
    
    vectorOutput[is.na(vectorOutput)] <-"NA"
   
    
    return(vectorOutput)
}

#-------------------------------------------------------------------------------
vectorOutputMatrices <- function(phenTestResult,outputMessages=TRUE){
    stop_message <- ""
    if (phenTestResult$method =="FE"){
        levels <-length(rownames(phenTestResult$model.output$count_matrix_all))
        if (levels>10){
            stop_message <- "Error:\nToo many levels for dependent variable.\n"
        }
        else{
            vectorOutput <- c(phenTestResult$depVariable, 
                    colnames(phenTestResult$model.output$count_matrix_all)[1], 
                    colnames(phenTestResult$model.output$count_matrix_all)[2],
                    rownames(phenTestResult$model.output$count_matrix_all))
            
            add_levels <- vector()
            
            if (levels<10){            
                add_levels<-rep(NA,each=(10-levels))
                vectorOutput<-c(vectorOutput,add_levels)
            }
            
            add_levels2 <-rep(NA,each=20)
            
            all <- as.vector(t(phenTestResult$model.output$count_matrix_all))
            
            vectorOutput <- c(vectorOutput,all,add_levels,add_levels)
            
            if (!is.null(phenTestResult$model.output$male)){
                males <- as.vector(t(phenTestResult$model.output$count_matrix_male))
                vectorOutput <- c(vectorOutput,males,add_levels,add_levels)
            }
            else{
                vectorOutput <- c(vectorOutput,add_levels2)  
            }   
            if (!is.null(phenTestResult$model.output$female)){
                females <- as.vector(t(phenTestResult$model.output$count_matrix_female))
                vectorOutput <- c(vectorOutput,females,add_levels,add_levels)
            }
            else {
                vectorOutput <- c(vectorOutput,add_levels2)  
            }
            
            names(vectorOutput) <- c("Dependent variable",
                    "Gp1 Genotype (g1)",
                    "Gp2 Genotype (g2)",
                    "Dependent variable level1 (l1)",
                    "Dependent variable level2 (l2)",
                    "Dependent variable level3 (l3)",
                    "Dependent variable level4 (l4)",
                    "Dependent variable level5 (l5)",
                    "Dependent variable level6 (l6)",
                    "Dependent variable level7 (l7)",
                    "Dependent variable level8 (l8)",
                    "Dependent variable level9 (l9)",
                    "Dependent variable level10 (l10)",
                    "Value g1_l1","Value g2_l1",
                    "Value g1_l2","Value g2_l2",
                    "Value g1_l3","Value g2_l3",
                    "Value g1_l4","Value g2_l4",
                    "Value g1_l5","Value g2_l5",
                    "Value g1_l6","Value g2_l6",
                    "Value g1_l7","Value g2_l7",
                    "Value g1_l8","Value g2_l8",
                    "Value g1_l9","Value g2_l9",
                    "Value g1_l10","Value g2_l10",
                    "Male Value g1_l1",
                    "Male Value g2_l1",
                    "Male Value g1_l2",
                    "Male Value g2_l2",
                    "Male Value g1_l3",
                    "Male Value g2_l3",
                    "Male Value g1_l4",
                    "Male Value g2_l4",
                    "Male Value g1_l5",
                    "Male Value g2_l5",
                    "Male Value g1_l6",
                    "Male Value g2_l6",
                    "Male Value g1_l7",
                    "Male Value g2_l7",
                    "Male Value g1_l8",
                    "Male Value g2_l8",
                    "Male Value g1_l9",
                    "Male Value g2_l9",
                    "Male Value g1_l10",
                    "Male Value g2_l10",
                    "Female Value g1_l1",
                    "Female Value g2_l1",
                    "Female Value g1_l2",
                    "Female Value g2_l2",
                    "Female Value g1_l3",
                    "Female Value g2_l3",
                    "Female Value g1_l4",
                    "Female Value g2_l4",
                    "Female Value g1_l5",
                    "Female Value g2_l5",
                    "Female Value g1_l6",
                    "Female Value g2_l6",
                    "Female Value g1_l7",
                    "Female Value g2_l7",
                    "Female Value g1_l8",
                    "Female Value g2_l8",
                    "Female Value g1_l9",
                    "Female Value g2_l9",
                    "Female Value g1_l10",
                    "Female Value g2_l10")
            return (vectorOutput)
        }
    } 
    else {
        stop_message <- "Error:\nThere are no matrices to output for MM framework. Function returns result only for FE framework.\n"
    }   
    
    if (nchar(stop_message)>0){
        if (outputMessages){
            message(stop_message)
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt))
            stop()
        }
        else {
            stop(stop_message)
        }
    } 
}
##------------------------------------------------------------------------------