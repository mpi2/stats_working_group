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
## JSONOutput.R contains JSONOutput function that returns analysis results in JSON format
##------------------------------------------------------------------------------
## Wrapper to prepare the output of the modeling and testing results in JSON 
## format. Assumes that modeling results are stored in the phenTestResult object 
## (output from functions testDataset and buildFinalModel)
JSONOutput <- function(phenTestResult, phenotypeThreshold=0.01)
{
    depVariable <- getVariable(phenTestResult)
    analysisResults <- analysisResults(phenTestResult)       
    noSexes <- length(levels(analysedDataset(phenTestResult)$Sex))
    
    refGenotypeVal <- paste('"reference genotype":"',refGenotype(phenTestResult),'"',sep="")
    testGenotypeVal <- paste('"test genotype":"',testGenotype(phenTestResult),'"',sep="")
    variable <- paste('"variable":"',depVariable,'"',sep="")
    
    if (method(phenTestResult) %in% c("MM","TF")) {        
        equation <- switch(analysisResults$equation,
                withoutWeight = {'"equation": "equation without weight"'},withWeight = {'"equation":"equation with weight"'})
        
        framework <- switch(method(phenTestResult),
                MM = '"method":"Mixed Model framework"',
                TF = '"method":"Time as Fixed Effect framework"')
        
        fittingMethod <- '"fittingMethod":"generalized least squares"'
        if (analysisResults$model.effect.batch && method(phenTestResult)=="MM")
        fittingMethod <- '"fittingMethod":"linear mixed-effects model"'        
        
        classificationValue <- paste('"classification tag":','"',classificationTag(phenTestResult,
                userMode="vectorOutput",outputMessages=FALSE),'"',sep="")

        
        x <- analysedDataset(phenTestResult)
        columnOfInterest <- x[,c(depVariable)]
        variability = paste('"variability of ',depVariable,' values":', 
                round(length(unique(columnOfInterest))/length(columnOfInterest),digits=3),sep="")
        
        Genotype_levels=levels(x$Genotype)
        Sex_levels=levels(x$Sex)
        
        DSsize = ""
        
        for (i in 1:length(Genotype_levels)){
            GenotypeSubset <- subset(x, x$Genotype==Genotype_levels[i])
            for (j in 1:length(Sex_levels)){
                GenotypeSexSubset <- subset(GenotypeSubset, GenotypeSubset$Sex==Sex_levels[j])
                columnOfInterest <- GenotypeSexSubset[,c(depVariable)]
                DSsize = paste(DSsize,paste('"',paste(Sex_levels[j],"(",Genotype_levels[i],")",sep=""),'"',sep=""),sep="")
                DSsize = paste(DSsize,":",sep="")
                DSsize = paste(DSsize,length(columnOfInterest),",",sep="")
            }
        }
        
        formula <- paste('"formula":','"',paste(format(analysisResults$model.formula.genotype), collapse= ' '),'"',sep="")
        
        
        
        addInfo = paste('"dataset details":{',DSsize,variability,"}",sep="")
        
        percentageChanges <- "NA"
        if (noSexes==2){
            percentageChanges <- paste('"Female":"',round(analysisResults$model.output.percentageChanges[1],digits=2),'%"',
                    ", ",
                    '"Male":"',
                    round(analysisResults$model.output.percentageChanges[2],digits=2),'%"',
                    sep="")
        } 
        else {
            if ("Female" %in% levels(x$Sex)){
                percentageChanges <- paste('"Female":"',round(analysisResults$model.output.percentageChanges[1],digits=2),'%"',
                        ", ",
                        '"Male":"NA"',
                        sep="")
                
            }
            else{
                percentageChanges <- paste('"Female":"NA"',
                        ", ",
                        '"Male":"',round(analysisResults$model.output.percentageChanges[1],digits=2),'%"',
                        sep="")
                
            }
            
        }   
        percentageChanges <- paste('"genotype percentage change":{',percentageChanges,"}",sep="") 
        transformation <- paste('"transformation":{',transformationJSON(phenTestResult),"}",sep="")
        genotype_pValMain <- paste('"p-value from the test of genotype contribution":',analysisResults$model.output.genotype.nulltest.pVal,sep="")
        
        statisticalAnalysis <- paste('"statistical analysis": {',
                paste(framework,fittingMethod,
                refGenotypeVal,testGenotypeVal,
                #equation,
                variable,
                transformation,
                formula,
                genotype_pValMain,
                percentageChanges,
                classificationValue,sep=","),"}",sep="")
        
        batch <- paste('"batch included":"',analysisResults$model.effect.batch,'"',sep="")
        weight <- paste('"weight included":"',analysisResults$model.effect.weight,'"',sep="")
        interaction <- paste('"sex and genotype interaction included":"',analysisResults$model.effect.interaction,'"',sep="")
        interaction_pval <- paste('"sex and genotype interaction effect p-value":',
                analysisResults$model.output.interaction,sep="")
        variance <- paste('"residual variances homogeneity":"',analysisResults$model.effect.variance,'"',sep="")

        
        genotype_estimate<- paste('"estimated coefficient":',
                analysisResults$model.output.summary["genotype_estimate"],sep="")
        genotype_se <- paste('"standard error":',
                analysisResults$model.output.summary["genotype_estimate_SE"],sep="")
        genotype_pVal <- paste('"p-value (represents the probability that effect is NOT relevant)":',
                analysisResults$model.output.summary["genotype_p_value"],sep="")
        sex_estimate<- paste('"estimated coefficient":',
                analysisResults$model.output.summary["sex_estimate"],sep="")
        sex_se <- paste('"standard error":',
                analysisResults$model.output.summary["sex_estimate_SE"],sep="")
        sex_pVal <- paste('"p-value (represents the probability that effect is NOT relevant)":',
                analysisResults$model.output.summary["sex_p_value"],sep="")
        weight_estimate<- paste('"estimated coefficient":',
                analysisResults$model.output.summary["weight_estimate"],sep="")
        weight_se <- paste('"standard error":',
                analysisResults$model.output.summary["weight_estimate_SE"],sep="")
        weight_pVal <- paste('"p-value (represents the probability that effect is NOT relevant)":',
                analysisResults$model.output.summary["weight_p_value"],sep="")
        
        intercept_estimate<- paste('"estimated coefficient":',
                analysisResults$model.output.summary["intercept_estimate"],sep="")
        intercept_se <- paste('"standard error":',
                analysisResults$model.output.summary["intercept_estimate_SE"],sep="")
        
        sexFvKO_estimate<- paste('"estimated coefficient":',
                analysisResults$model.output.summary["sex_FvKO_estimate"],sep="")
        sexFvKO_se <- paste('"standard error":',
                analysisResults$model.output.summary["sex_FvKO_SE"],sep="")
        sexFvKO_pVal <- paste('"p-value (represents the probability that effect is NOT relevant)":',
                analysisResults$model.output.summary["sex_FvKO_p_value"],sep="")
        sexMvKO_estimate<- paste('"estimated coefficient":',
                analysisResults$model.output.summary["sex_MvKO_estimate"],sep="")
        sexMvKO_se <- paste('"standard error":',
                analysisResults$model.output.summary["sex_MvKO_SE"],sep="")
        sexMvKO_pVal <- paste('"p-value (represents the probability that effect is NOT relevant)":',
                analysisResults$model.output.summary["sex_MvKO_p_value"],sep="")
        
        effects <- paste('"effects": {',
                '"intercept":{',paste(intercept_estimate,intercept_se,sep=","),"},",
                '"genotype":{',paste(genotype_estimate,genotype_se,genotype_pVal,sep=","),"},",
                '"sex":{',paste(sex_estimate,sex_se,sex_pVal,sep=","),"},",
                '"weight":{',paste(weight_estimate,weight_se,weight_pVal,sep=","),"},",
                '"females in test genotype subset":{',paste(sexFvKO_estimate,sexFvKO_se,sexFvKO_pVal,sep=","),"},",
                '"males in test genotype subset":{',paste(sexMvKO_estimate,sexMvKO_se,sexMvKO_pVal,sep=","),"}",
                "}",sep="")
        


        quality <- paste('"residuals normality test": {"reference genotype group":',
                analysisResults$model.output.quality[[2]],',',
                '"test genotype group":',
                analysisResults$model.output.quality[[4]],'},',
                '"blups test":',
                analysisResults$model.output.quality[[5]],',',
                '"rotated residuals normality test":',
                analysisResults$model.output.quality[[6]],sep="")
        
        quality <- paste('"model fitting quality": {',quality,'}',sep="")

        statisticalAnalysisDetails <- paste('"details": {',
                paste(batch,weight,interaction,variance,effects,quality
                        ,sep=","),"}",sep="")
        
        JSONOutput <- paste(statisticalAnalysis,
                statisticalAnalysisDetails,
                addInfo,sep=",")
        
        
    }
    else if (method(phenTestResult) %in% c("FE","RR")){
        male_pval <- NA
        female_pval <- NA
        male_ES <- NA
        female_ES <- NA
        male <- NA
        female <- NA
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

            if (noSexes==2){                
                female <- paste('"females only dataset": {"high vs normal/low":{"p_value":',high_female_p.value,',"effect size":"',high_female_ES,'"},',
                                '"low vs normal/high":{"p_value":',low_female_p.value,',"effect size":"',low_female_ES,'"}}',sep="")
                male <- paste('"males only dataset": {"high vs normal/low":{"p_value":',high_male_p.value,',"effect size":"',high_male_ES,'"},',
                              '"low vs normal/high":{"p_value":',low_male_p.value,',"effect size":"',low_male_ES,'"}}',sep="")
            }
           
            
            all <- paste('"combined dataset": {"high vs normal/low":{"p_value":',high_all_p.value,',"effect size":"',high_all_ES,'"},',
                         '"low vs normal/high":{"p_value":',low_all_p.value,',"effect size":"',low_all_ES,'"}}',sep="")

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
            
            all <- paste('"combined dataset":{"p_value":',p_value_all,',"effect size":"',ES_all,'"}',sep="")
            female <- paste('"females only dataset":{"p_value":',female_pval,',"effect size":"',female_ES,'"}',sep="")
            male <- paste('"males only dataset":{"p_value":',male_pval,',"effect size":"',male_ES,'"}',sep="")
        }
        
        
        classificationValue <- paste('"classification tag":"',classificationTag(phenTestResult,phenotypeThreshold=phenotypeThreshold,
                outputMessages=FALSE),'"',sep="")
        
        if (method(phenTestResult)=="RR"){
            addInfo = paste('"',gsub(":","",rownames(parameters(phenTestResult))[1]),'":"',
                    parameters(phenTestResult)[1],'"',",",sep="")
            addInfo = paste(addInfo,'"',gsub(":","",rownames(parameters(phenTestResult))[2]),'":"',
                    parameters(phenTestResult)[2],'"',",",sep="")
            
            if (noSexes==2){
                addInfo = paste(addInfo,'"',gsub(":","",rownames(parameters(phenTestResult))[3]),'":"',
                        parameters(phenTestResult)[3],'"',",",sep="")
                addInfo = paste(addInfo,'"',gsub(":","",rownames(parameters(phenTestResult))[4]),'":"',
                        parameters(phenTestResult)[4],'"',sep="")
            }
            else {
                addInfo = paste(addInfo,'"',gsub(":","",rownames(parameters(phenTestResult))[3]),'":"',
                        parameters(phenTestResult)[3],'"',sep="")
            }
            classificationValue = paste(classificationValue,', "RR method details":{',addInfo,'}',sep="")
        }
        
        
        framework <- switch(method(phenTestResult),FE = '"method":"Fisher Exact Test framework"',
                RR = '"method":"Reference Ranges Plus framework"')

        if (!is.na(female)) all <- paste(all,female,male,sep=",")
        JSONOutput <- paste('"statistical analysis": {',
                paste(framework,refGenotypeVal,testGenotypeVal,
                        variable,all,
                        classificationValue,sep=","),"}",sep="")
        
        

    }
    else if (method(phenTestResult) == "LR"){        
        
        classificationValue <- paste('"classification tag":"',
                classificationTag(phenTestResult, userMode="vectorOutput",
                        outputMessages=FALSE),'"',sep="")
        
        framework <- '"method":"Logistic Regression framework"'
        genotype_pValMain <- paste('"genotype p-value":',analysisResults$model.output.genotype.nulltest.pVal,sep="")
        formula <- paste('"formula":','"',paste(format(analysisResults$model.formula.genotype), collapse= ' '),'"',sep="")
        
        statisticalAnalysis <- paste('"statistical analysis": {',
                paste(framework,formula,refGenotypeVal,testGenotypeVal,
                        variable,genotype_pValMain,
                        classificationValue,sep=","),"}",sep="")
        
        batch <- paste('"batch included":"',analysisResults$model.effect.batch,'"',sep="")
        interaction <- paste('"sex and genotype interaction included":"',analysisResults$model.effect.interaction,'"',sep="")
        interaction_pval <- paste('"sex and genotype interaction effect p-value":',
                analysisResults$model.output.interaction,sep="")
        
        genotype_estimate<- paste('"estimated coefficient":',
                analysisResults$model.output.summary["genotype_estimate"],sep="")
        genotype_se <- paste('"standard error":',
                analysisResults$model.output.summary["genotype_estimate_SE"],sep="")
        genotype_pVal <- paste('"p-value (represents the probability that effect is NOT relevant)":',
                analysisResults$model.output.summary["genotype_p_value"],sep="")
        sex_estimate<- paste('"estimated coefficient":',
                analysisResults$model.output.summary["sex_estimate"],sep="")
        sex_se <- paste('"standard error":',
                analysisResults$model.output.summary["sex_estimate_SE"],sep="")
        sex_pVal <- paste('"p-value (represents the probability that effect is NOT relevant)":',
                analysisResults$model.output.summary["sex_p_value"],sep="")
        
        intercept_estimate<- paste('"estimated coefficient":',
                analysisResults$model.output.summary["intercept_estimate"],sep="")
        intercept_se <- paste('"standard error":',
                analysisResults$model.output.summary["intercept_estimate_SE"],sep="")
        
        sexFvKO_estimate<- paste('"estimated coefficient":',
                analysisResults$model.output.summary["sex_FvKO_estimate"],sep="")
        sexFvKO_se <- paste('"standard error":',
                analysisResults$model.output.summary["sex_FvKO_SE"],sep="")
        sexFvKO_pVal <- paste('"p-value (represents the probability that effect is NOT relevant)":',
                analysisResults$model.output.summary["sex_FvKO_p_value"],sep="")
        sexMvKO_estimate<- paste('"estimated coefficient":',
                analysisResults$model.output.summary["sex_MvKO_estimate"],sep="")
        sexMvKO_se <- paste('"standard error":',
                analysisResults$model.output.summary["sex_MvKO_SE"],sep="")
        sexMvKO_pVal <- paste('"p-value (represents the probability that effect is NOT relevant)":',
                analysisResults$model.output.summary["sex_MvKO_p_value"],sep="")
        
        effects <- paste('"effects": {',
                '"intercept":{',paste(intercept_estimate,intercept_se,sep=","),"},",
                '"genotype":{',paste(genotype_estimate,genotype_se,genotype_pVal,sep=","),"},",
                '"sex":{',paste(sex_estimate,sex_se,sex_pVal,sep=","),"},",
                '"females in test genotype subset":{',paste(sexFvKO_estimate,sexFvKO_se,sexFvKO_pVal,sep=","),"},",
                '"males in test genotype subset":{',paste(sexMvKO_estimate,sexMvKO_se,sexMvKO_pVal,sep=","),"}",
                "}",sep="")
        
        
   
        
        statisticalAnalysisDetails <- paste('"details": {',
                paste(batch,interaction,effects
                        ,sep=","),"}",sep="")
        
        JSONOutput <- paste(statisticalAnalysis,
                statisticalAnalysisDetails,sep=",")
        

    }
    
    JSONOutput <- gsub("NA", "null", JSONOutput)
    JSONOutput <- gsub('"null"', "null", JSONOutput)
    return(paste("{",JSONOutput,"}",sep=""))
}

#-------------------------------------------------------------------------------
