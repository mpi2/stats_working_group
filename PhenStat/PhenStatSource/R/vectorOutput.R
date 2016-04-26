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
## vectorOutput.R contains vectorOutput and outputLength functions
##------------------------------------------------------------------------------
## Wrapper to prepare the output of the modeling and testing results in vector 
## form. Assumes that modeling results are stored in the phenTestResult object 
## (output from functions testDataset and buildFinalModel)
vectorOutput <- function(phenTestResult, phenotypeThreshold=0.01)
{
    depVariable <- getVariable(phenTestResult)
    analysisResults <- analysisResults(phenTestResult)       
    #noSexes <- length(levels(analysedDataset(phenTestResult)$Sex))
    
    if (method(phenTestResult) %in% c("MM","TF")) {        
        equation <- switch(analysisResults$equation,
                withoutWeight = {"equation withoutWeight"},withWeight = {"equation withWeight"})
        
        framework <- switch(method(phenTestResult),
                MM = "Mixed Model framework",
                TF = "Time as Fixed Effect framework")
        
        fittingMethod <- "generalized least squares, "
        if (analysisResults$model.effect.batch && method(phenTestResult)=="MM")
        fittingMethod <- "linear mixed-effects model, "        
        
        classificationValue <- classificationTag(phenTestResult,
                userMode="vectorOutput",outputMessages=FALSE)
        
        
        x <- analysedDataset(phenTestResult)
        columnOfInterest <- x[,c(depVariable)]
        variability = paste('"variability":', round(length(unique(columnOfInterest))/length(columnOfInterest),digits=3),sep="")
        
        Genotype_levels=levels(x$Genotype)
        Sex_levels=levels(x$Sex)
        
        DSsize = ""
        
        for (i in 1:length(Genotype_levels)){
            GenotypeSubset <- subset(x, x$Genotype==Genotype_levels[i])
            for (j in 1:length(Sex_levels)){
                GenotypeSexSubset <- subset(GenotypeSubset, GenotypeSubset$Sex==Sex_levels[j])
                columnOfInterest <- GenotypeSexSubset[,c(depVariable)]
                DSsize = paste(DSsize,paste('"',paste(Genotype_levels[i],Sex_levels[j],sep="_"),'"',sep=""),sep="")
                DSsize = paste(DSsize,":",sep="")
                DSsize = paste(DSsize,length(columnOfInterest),",",sep="")
            }
        }
        
        formula <- paste('"Formula":','"',paste(format(analysisResults$model.formula.genotype), collapse= ' '),'"',sep="")
        
        
        
        addInfo = paste("{",DSsize,variability,",",formula,"}",sep="")
        
        percentageChanges <- "NA"
        if (noSexes(phenTestResult)==2){
            percentageChanges <- paste("Female: ",round(analysisResults$model.output.percentageChanges[1],digits=2),"%",
                    ", ",
                    "Male: ",
                    round(analysisResults$model.output.percentageChanges[2],digits=2),"%",
                    sep="")
        } 
        else {
            if ("Female" %in% levels(x$Sex)){
                percentageChanges <- paste("Female: ",round(analysisResults$model.output.percentageChanges[1],digits=2),"%",
                        ", ",
                        "Male: NA",
                        sep="")
                
            }
            else{
                percentageChanges <- paste("Female: NA",
                        ", ",
                        "Male: ",round(analysisResults$model.output.percentageChanges[1],digits=2),"%",
                        sep="")
                
            }
            
        }   
        
        vectorOutput <- c(paste(framework,", ",fittingMethod, equation,sep=""),
                as.character(depVariable), 
                as.character(analysisResults$model.effect.batch), 
                as.character(analysisResults$model.effect.variance), 
                as.character(analysisResults$model.output.genotype.nulltest.pVal), 
                as.character(analysisResults$model.output.summary["genotype_estimate"]), 
                as.character(analysisResults$model.output.summary["genotype_estimate_SE"]),  
                as.character(analysisResults$model.output.summary["genotype_p_value"]),
                as.character(percentageChanges),
                as.character(analysisResults$model.output.summary["sex_estimate"]), 
                as.character(analysisResults$model.output.summary["sex_estimate_SE"]),  
                as.character(analysisResults$model.output.summary["sex_p_value"]), 
                as.character(analysisResults$model.output.summary["weight_estimate"]), 
                as.character(analysisResults$model.output.summary["weight_estimate_SE"]), 
                as.character(analysisResults$model.output.summary["weight_p_value"]), 
                as.character(analysisResults$model.output.quality), 
                as.character(analysisResults$model.output.summary["intercept_estimate"]), 
                as.character(analysisResults$model.output.summary["intercept_estimate_SE"]), 
                as.character(analysisResults$model.effect.interaction),
                as.character(analysisResults$model.output.interaction),
                as.character(analysisResults$model.output.summary["sex_FvKO_estimate"]), 
                as.character(analysisResults$model.output.summary["sex_FvKO_SE"]), 
                as.character(analysisResults$model.output.summary["sex_FvKO_p_value"]),  
                as.character(analysisResults$model.output.summary["sex_MvKO_estimate"]),  
                as.character(analysisResults$model.output.summary["sex_MvKO_SE"]), 
                as.character(analysisResults$model.output.summary["sex_MvKO_p_value"]),
                as.character(classificationValue),
                transformation(phenTestResult),
                as.character(addInfo))
        
        names(vectorOutput) <- c("Method",
                "Dependent variable",
                "Batch included",
                "Residual variances homogeneity",
                "Genotype contribution",
                "Genotype estimate",
                "Genotype standard error",
                "Genotype p-Val",
                "Genotype percentage change",
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
                "Transformation",
                "Additional information")
    }
    else if (method(phenTestResult) %in% c("FE","RR")){
        male_pval <- NA
        female_pval <- NA
        male_ES <- NA
        female_ES <- NA
        if (method(phenTestResult)=="RR"){
            for (i in seq_along(analysisResults)) {
                val <- analysisResults[[i]]
                if (comparison(val)=="High vs Normal/Low") {
                    if (analysedSubset(val)=="all"){
                        high_all_p.value<-as.numeric(as.vector(getColumnView(val))[2]) 
                        high_all_ES<-as.character(as.vector(getColumnView(val))[3])       
                    }
                    if (analysedSubset(val)=="males"){
                        high_male_p.value<-as.numeric(as.vector(getColumnView(val))[2])  
                        high_male_ES<-as.character(as.vector(getColumnView(val))[3])                       
                    }
                    if (analysedSubset(val)=="females"){
                        high_female_p.value<-as.numeric(as.vector(getColumnView(val))[2])  
                        high_female_ES<-as.character(as.vector(getColumnView(val))[3])                       
                    }
                }
                else {
                    if (analysedSubset(val)=="all"){
                        low_all_p.value<-as.numeric(as.vector(getColumnView(val))[2])     
                        low_all_ES<-as.character(as.vector(getColumnView(val))[3])                    
                    }
                    if (analysedSubset(val)=="males"){
                        low_male_p.value<-as.numeric(as.vector(getColumnView(val))[2])  
                        low_male_ES<-as.character(as.vector(getColumnView(val))[3])                    
                    }
                    if (analysedSubset(val)=="females"){
                        low_female_p.value<-as.numeric(as.vector(getColumnView(val))[2])       
                        low_female_ES<-as.character(as.vector(getColumnView(val))[3])               
                    }
                }
            }
            p_value_all <- paste(low_all_p.value,high_all_p.value,sep=",")
            ES_all <- paste(low_all_ES,high_all_ES,sep=",")
            if (noSexes(phenTestResult)==2){
                male_pval <- paste(low_male_p.value,high_male_p.value,sep=",")
                female_pval <- paste(low_female_p.value,high_female_p.value,sep=",")
                male_ES <- paste(low_male_ES,high_male_ES,sep=",")
                female_ES <- paste(low_female_ES,high_female_ES,sep=",")
            }
        }
        if (method(phenTestResult)=="FE"){
            for (i in seq_along(analysisResults)) {
                val <- analysisResults[[i]]
                if (analysedSubset(val)=="all"){
                    p_value_all<-as.numeric(as.vector(getColumnView(val))[2]) 
                    ES_all<-as.character(as.vector(getColumnView(val))[3])       
                }
                if (analysedSubset(val)=="males"){
                    male_pval<-as.numeric(as.vector(getColumnView(val))[2])  
                    male_ES<-as.character(as.vector(getColumnView(val))[3])                       
                }
                if (analysedSubset(val)=="females"){
                    female_pval<-as.numeric(as.vector(getColumnView(val))[2])  
                    female_ES<-as.character(as.vector(getColumnView(val))[3])                       
                }
            }
        }
        
      
        classificationValue <- classificationTag(phenTestResult,phenotypeThreshold=phenotypeThreshold,
                outputMessages=FALSE)
        
        if (method(phenTestResult)=="RR"){
            addInfo = paste('"',rownames(parameters(phenTestResult))[1],
                    parameters(phenTestResult)[1],'"',",",sep="")
            addInfo = paste(addInfo,'"',rownames(parameters(phenTestResult))[2],
                            parameters(phenTestResult)[2],'"',",",sep="")
            
            if (noSexes(phenTestResult)==2){
                addInfo = paste(addInfo,'"',rownames(parameters(phenTestResult))[3],
                                parameters(phenTestResult)[3],'"',",",sep="")
                addInfo = paste(addInfo,'"',rownames(parameters(phenTestResult))[4],
                                parameters(phenTestResult)[4],'"',sep="")
            }
            else {
                addInfo = paste(addInfo,'"',rownames(parameters(phenTestResult))[3],
                                parameters(phenTestResult)[3],'"',sep="")
            }
            addInfo = paste("{",addInfo,"}",sep="")
        }
        else {
            addInfo = "NA"
        }
        
        
        vectorOutput <- c(switch(method(phenTestResult),FE = "Fisher Exact Test framework",
                        RR = "Reference Ranges Plus framework"),
                as.character(depVariable), 
                "NA", 
                "NA",
                "NA",
                ES_all,
                "NA",  
                p_value_all,
                "NA",
                "NA",
                "NA", #10 
                "NA", 
                "NA", 
                "NA", 
                "NA", 
                as.character(testGenotype(phenTestResult)), 
                "NA",#16 
                as.character(refGenotype(phenTestResult)),  
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
                "lambda=NA, scaleShift=NA",
                addInfo)
        
        names(vectorOutput) <- c("Method",
                "Dependent variable",
                "Batch included",
                "Residual variances homogeneity",
                "Genotype contribution",
                "Genotype estimate",
                "Genotype standard error",
                "Genotype p-Val",
                "Genotype percentage change",
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
                "Transformation",
                "Additional information")
    }
    else if (method(phenTestResult) == "LR"){        
        
        classificationValue <- classificationTag(phenTestResult, userMode="vectorOutput",outputMessages=FALSE)
        
        vectorOutput <- c(
                "Logistic Regression framework",
                as.character(depVariable), 
                as.character(analysisResults$model.effect.batch), 
                "NA",  
                as.character(analysisResults$model.output.genotype.nulltest.pVal), 
                as.character(analysisResults$model.output.summary["genotype_estimate"]), 
                as.character(analysisResults$model.output.summary["genotype_estimate_SE"]),  
                as.character(analysisResults$model.output.summary["genotype_p_value"]),
                "NA" , 
                as.character(analysisResults$model.output.summary["sex_estimate"]), 
                as.character(analysisResults$model.output.summary["sex_estimate_SE"]),  
                as.character(analysisResults$model.output.summary["sex_p_value"]), 
                "NA", 
                "NA", 
                "NA",
                as.character(analysisResults$model.output.quality), 
                as.character(analysisResults$model.output.summary["intercept_estimate"]), 
                as.character(analysisResults$model.output.summary["intercept_estimate_SE"]), 
                as.character(analysisResults$model.effect.interaction),
                as.character(analysisResults$model.output.interaction),
                as.character(analysisResults$model.output.summary["sex_FvKO_estimate"]), 
                as.character(analysisResults$model.output.summary["sex_FvKO_SE"]), 
                as.character(analysisResults$model.output.summary["sex_FvKO_p_value"]),  
                as.character(analysisResults$model.output.summary["sex_MvKO_estimate"]),  
                as.character(analysisResults$model.output.summary["sex_MvKO_SE"]), 
                as.character(analysisResults$model.output.summary["sex_MvKO_p_value"]),
                as.character(classificationValue),
                "lambda=NA, scaleShift=NA",
                "NA" )  
        
        names(vectorOutput) <- c("Method",
                "Dependent variable",
                "Batch included",
                "Residual variances homogeneity",
                "Genotype contribution",
                "Genotype estimate",
                "Genotype standard error",
                "Genotype p-Val",
                "Genotype percentage change",
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
                "Transformation",
                "Additional information")
    }
    
    vectorOutput[is.na(vectorOutput)] <-"NA"
    
    
    return(vectorOutput)
}

#-------------------------------------------------------------------------------
vectorOutputMatrices <- function(phenTestResult,outputMessages=TRUE){
    stop_message <- ""
    count_matrices <- getCountMatrices(phenTestResult)
    #noSexes <- length(levels(analysedDataset(phenTestResult)$Sex))
    if (method(phenTestResult) %in% c("FE","RR")){
        if (method(phenTestResult)=="FE"){
          values <- analysedDataset(phenTestResult)[, getVariable(phenTestResult)]
          values <- factor(values)
          levels <- length(levels(values))
        }
        else {
          levels <- length(rownames(count_matrices[["all"]]))
        }

        vectorOutput <- c(getVariable(phenTestResult), 
                colnames(count_matrices[["all"]])[1], 
                colnames(count_matrices[["all"]])[2],
                rownames(count_matrices[["all"]]))
        
        add_levels <- vector()
        
        if (levels<10){            
            add_levels<-rep(NA,each=(10-levels))
            vectorOutput<-c(vectorOutput,add_levels)
        }
        
        add_levels2 <-rep(NA,each=20)
        
        all <- as.vector(t(count_matrices[["all"]]))
        
        vectorOutput <- c(vectorOutput,all,add_levels,add_levels)
        
        if (noSexes(phenTestResult)==2){
            males <- as.vector(t(count_matrices[["male"]]))
            vectorOutput <- c(vectorOutput,males,add_levels,add_levels)
        }
        else{
            vectorOutput <- c(vectorOutput,add_levels2)  
        }   
        if (noSexes(phenTestResult)==2){
            females <- as.vector(t(count_matrices[["female"]]))
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
    else {
        stop_message <- paste("Error:\nThere are no matrices to output for MM, TF or LR frameworks.", 
    "Function returns result only for FE and RR frameworks.\n",sep="")
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