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
                withoutWeight = {"Eq1"},withWeight = {"Eq2"})
        
        classificationValue <- classificationTag(phenTestResult,
                userMode="vectorOutput",outputMessages=FALSE)
        
        vectorOutput <- c(paste("MM - ",equation,sep=""),
                phenTestResult$depVariable, 
                phenTestResult$model.effect.batch, 
                phenTestResult$model.effect.variance, 
                phenTestResult$model.output.genotype.nulltest.pVal, 
                phenTestResult$model.output.summary["genotype_estimate"], 
                phenTestResult$model.output.summary["genotype_estimate_SE"],  
                phenTestResult$model.output.summary["genotype_p_value"],
                phenTestResult$model.output.summary["gender_estimate"], 
                phenTestResult$model.output.summary["gender_estimate_SE"],  
                phenTestResult$model.output.summary["gender_p_value"], 
                phenTestResult$model.output.summary["weight_estimate"], 
                phenTestResult$model.output.summary["weight_estimate_SE"], 
                phenTestResult$model.output.summary["weight_p_value"], 
                phenTestResult$model.output.quality, 
                phenTestResult$model.output.summary["intercept_estimate"], 
                phenTestResult$model.output.summary["intercept_estimate_SE"], 
                phenTestResult$model.effect.interaction,
                phenTestResult$model.output.interaction,
                phenTestResult$model.output.summary["gender_FvKO_estimate"], 
                phenTestResult$model.output.summary["gender_FvKO_SE"], 
                phenTestResult$model.output.summary["gender_FvKO_p_value"],  
                phenTestResult$model.output.summary["gender_MvKO_estimate"],  
                phenTestResult$model.output.summary["gender_MvKO_SE"], 
                phenTestResult$model.output.summary["gender_MvKO_p_value"],
                classificationValue)
        names(vectorOutput) <- c("Method",
                "Dependent variable",
                "Batch included",
                "Residual variances homogeneity",
                "Genotype contribution",
                "Genotype estimate",
                "Genotype standard error",
                "Genotype p-Val",
                "Gender estimate",
                "Gender standard error",
                "Gender p-val",
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
                "Gender FvKO estimate",
                "Gender FvKO standard error",
                "Gender FvKO p-val",
                "Gender MvKO estimate",
                "Gender MvKO standard error",
                "Gender MvKO p-val",
                "Classification tag")
    }
    else if (phenTestResult$method=="FE"){
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
        
        vectorOutput <- c("Fisher Exact Test",
                phenTestResult$depVariable, 
                NA, 
                NA,
                NA,
                as.numeric(phenTestResult$model.output$ES),
                NA,  
                as.numeric(phenTestResult$model.output$all$p.val),
                NA, 
                NA, #10 
                NA, 
                NA, 
                NA, 
                NA, 
                colnames(phenTestResult$model.output$count_matrix_all)[1], 
                NA,#16 
                colnames(phenTestResult$model.output$count_matrix_all)[2],  
                NA, 
                NA,
                NA, 
                NA, 
                NA,
                NA,
                NA, #24
                female_ES,  
                NA,  
                female_pval, 
                male_ES,  
                NA,  
                male_pval, #30
                classificationTag(phenTestResult))
        
        names(vectorOutput) <- c("Method",
                "Dependent variable",
                "Batch included",
                "Residual variances homogeneity",
                "Genotype contribution",
                "Genotype estimate",
                "Genotype standard error",
                "Genotype p-Val",
                "Gender estimate",
                "Gender standard error",
                "Gender p-val", #11
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
                "Gender FvKO estimate",
                "Gender FvKO standard error",
                "Gender FvKO p-val",
                "Gender MvKO estimate",
                "Gender MvKO standard error",
                "Gender MvKO p-val",
                "Classification tag")
    }
    return(vectorOutput)
}

#-------------------------------------------------------------------------------
vectorOutputMatrices <- function(phenTestResult){
    if (phenTestResult$method=="FE"){
        levels <-length(rownames(phenTestResult$model.output$count_matrix_all))
        if (levels>10){
            stop("Too many levels for dependent variable")
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
}
##------------------------------------------------------------------------------