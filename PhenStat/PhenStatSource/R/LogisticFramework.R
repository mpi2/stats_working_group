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
## LogisticFramework.R contains LRDataset, startLRModel, finalLRModel, modelFormulaLR
## and parserOutputSummaryLR functions
##------------------------------------------------------------------------------
## LRDataset prepares dataset for the LR framework - maps values to 0/1
## Abnormal - 1 (will be modeled), Normal - 0 
LRDataset <- function(phenList=NULL, depVariable=NULL, abnormalValues=c("abnormal","Abnormal","TRUE","deviant"),
        outputMessages=TRUE)
{
    stop_message <- ""
    ## START CHECK ARGUMENTS   
    # 1
    if (is.null(phenList) || !is(phenList,"PhenList")){
        stop_message <- paste(stop_message,"Error:\nPlease create and specify PhenList object.\n",sep="")
    }
    # 2
    if (is.null(depVariable)){
        stop_message <- paste(stop_message,"Error:\nPlease specify dependent variable, ",
                "for example: depVariable='Lean.Mass'.\n",sep="")
    }
    if (nchar(stop_message)==0) {
        # create small data frame with values to analyse
        columnOfInterest <- getColumn(phenList,depVariable)
        # Presence
        if (is.null(columnOfInterest))
            stop_message <- paste("Error:\nDependent variable column '",
                depVariable,"' is missed in the dataset.\n",sep="")
    }    
    ## STOP CHECK ARGUMENTS   
    # If problems are deteckted
    if (nchar(stop_message)>0){
        if (outputMessages)  { 
            message(stop_message)
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt))      
            stop()
        }
        else {
            stop(stop_message)
        }
    }
    
    # Creates a subset with analysable variable
    if (nchar(stop_message)==0) {
        oldvalues <- levels(factor(na.omit(columnOfInterest)))
        oldvaluesAbnormal <- oldvalues[which(oldvalues %in% abnormalValues)]
        oldvaluesNormal <- oldvalues[which(!(oldvalues %in% abnormalValues))]
        newvalues <- c(rep(1,length(oldvaluesAbnormal)),rep(0,length(oldvaluesNormal)))
        map <- setNames(newvalues,c(oldvaluesAbnormal,oldvaluesNormal))
        mappedvalues <- sapply(columnOfInterest, function(x) map[as.character(x)])
        names(mappedvalues)<-NULL
        columnOfInterest <- as.integer(mappedvalues)
        phenList@datasetPL[,c(depVariable)] <- columnOfInterest
        return(phenList)
    }  
}
##------------------------------------------------------------------------------
## LRstartModel creates start model and modifies it after testing of different
## hypothesis (the model effects).
## The model effects are:
## - batch effect (random effect significance)
## Fixed effects:
## - interaction effect (genotype by sex interaction significance)
## - sex effect (sex significance)
##
## As the result PhenTestResult object that contains calculated or user
## defined model effects and MM start model is created.
##http://idiom.ucsd.edu/~rlevy/lign251/fall2007/lecture_15.pdf 
## States for nested mixed effect models liklihood ratio can be used to compare model
##provided only differ in fixed effects structure
startLRModel <- function(phenList, depVariable, outputMessages=TRUE, pThreshold=0.05)
##------------------------------------------------------------------------------
{
    x <- getDataset(phenList)
    numberofsexes <- length(levels(x$Sex))
    keep_batch <- FALSE
    keep_interaction <- FALSE
    keep_sex <- FALSE
    keep_equalvar <- NA  
    keep_weight <- NA  # weight not included as a fixed effect      
    equation <- "withoutWeight"   # weight not included as a fixed effect    
    #Test for batch effect
    formula_withBatch=modelFormulaLR(numberofsexes, depVariable, sexIncluded=TRUE, dimorphismIncluded=TRUE, 
            IncludeBatch="Yes")
    formula_withOutBatch=modelFormulaLR(numberofsexes, depVariable, sexIncluded=TRUE, dimorphismIncluded=TRUE, 
            IncludeBatch="No")        
    suppressWarnings(
        try(## Goal of this section is to assess whether batch is significant or not in explaining variation
                    {
        if (batchIn(phenList)){
                            ## LR fit of model formula with no random effects
                            ## Model 1A (model_withoutBatch)
                            ## glm - Fit a generalized linear models, can not handle nested effects.
                            L_model_withoutbatch <- do.call("glm",args=list(formula_withOutBatch, x, 
                                            na.action="na.omit", family=binomial()  ))
                            ## LR fit of model formula (with random effects)
                            ## glmer - Fit a generalized linear mixed model (GLMM) can linear mixed effect models 
                            ## with fixed effects and random/nested effects.                            
                            ## Model 1 (model_withBatch)
                            L_model_withBatch <- do.call("glmer", args = list(formula_withBatch, data = x, 
                                            na.action="na.omit",family=binomial() ))    
                            ## Test: the random effects associated with batch intercepts can be ommited from model
                            ## Hypothesis 1
                            ## Null Hypothesis: variance of batch = 0
                            ## Alternative Hypothesis: variance of batch > 0    
                            ## If p value below threshold then you reject null and accept alternative that batch is 
                            ## significant in explaining variation in the model
                            ## Based on method shown here: 
                            ## http://www.simonqueenborough.com/R/specialist/mixed-models.html
                            p.value.batch <- pchisq(-2*(logLik(L_model_withBatch)-logLik(L_model_withoutbatch)), 
                                    1, lower.tail=FALSE)[1] 
                            keep_batch <- p.value.batch<pThreshold
                            
                        }                                             
                    },
                   silent=TRUE)
    )
    
    #START OF tryCatch    
    finalResult <- tryCatch({                
                ## Goal of this section is to tests for significance of fixed effects by comparing models 
                ## with and without fixed effects included        
                if(numberofsexes==2){                    
                    #test sexual dimorphism
                    L_model_withoutbatch <- do.call("logistf",args=list(formula_withOutBatch, x, na.action="na.omit"))
                    
                    formula_withoutBatch_noSD=modelFormulaLR(numberofsexes, depVariable, sexIncluded=TRUE, 
                            dimorphismIncluded=FALSE, IncludeBatch="No")
                    L_model_withoutBatch_noSD <- do.call("logistf", args = list(formula_withoutBatch_noSD, data = x, 
                                    na.action="na.omit" ))                    
                    interactionTest = anova(L_model_withoutbatch, L_model_withoutBatch_noSD)$pval
                    keep_interaction <- interactionTest<pThreshold
                                        
                    #test sex inclusion in model
                    formula_withoutBatch_noSex=modelFormulaLR(numberofsexes, depVariable, sexIncluded=FALSE, 
                            dimorphismIncluded=FALSE, IncludeBatch="No")
                    L_model_withoutBatch_noSex <- do.call("logistf", args = list(formula_withoutBatch_noSex, 
                                    data = x, na.action="na.omit"))                    
                    sexTest = anova(L_model_withoutbatch, L_model_withoutBatch_noSex)$pval
                    keep_sex <- sexTest<pThreshold
                }
                else {
                    keep_sex <- FALSE
                    keep_interaction <- FALSE
                    interactionTest <- NA
                }
                
                if (outputMessages)
                message(paste("Information:\nCalculated values for model effects are: ",
                                "keepBatch=",keep_batch,
                                ", keepSex=",keep_sex,
                                ", keepInteraction=",keep_interaction,".\n",sep=""))                
                ## Need to obtain final model 
                ## step one build formula and then fit the model through glm or glmer route depending on 
                ## whether batch is significant          
                finalFormula=modelFormulaLR(numberofsexes, depVariable, sexIncluded=keep_sex, 
                        dimorphismIncluded=keep_interaction, IncludeBatch="No")    
                finalModel <- do.call("logistf", args = list(finalFormula, data = x, na.action="na.omit"))
                #LR modelling output 
                linearRegressionOutput <- list(model.output=finalModel,
                        equation = equation,
                        model.effect.batch=keep_batch,
                        model.effect.interaction=keep_interaction,
                        model.output.interaction=interactionTest,
                        model.effect.sex=keep_sex,
                        model.effect.weight=keep_weight,
                        model.effect.variance=keep_equalvar,
                        numberSexes=numberofsexes,
                        pThreshold=pThreshold,
                        model.formula.genotype=finalFormula)
                linearRegressionOutput$FET <- FisherExactTest(phenList,depVariable) 
                
                 finalResult <- new("PhenTestResult",
                    analysedDataset=x,
                    depVariable=depVariable,
                    refGenotype=refGenotype(phenList),
                    testGenotype=testGenotype(phenList),
                    method="LR",
                    transformationRequired = FALSE,
                    lambdaValue = integer(0),
                    scaleShift = integer(0),
                    analysisResults=linearRegressionOutput) 
            },
            #END OF tryCatch    
            error=function(error_mes) {
                finalResult <- NULL                
                stop_message <- paste("Error:\nCan't fit the model ",
                        format(formula_withOutBatch),". Try FE method.\n",sep="")                  
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
    )                    
    return(finalResult)    
}
##------------------------------------------------------------------------------
## Function to build start formula for Logistic framework 
modelFormulaLR <- function(numberofsexes, depVariable, sexIncluded, dimorphismIncluded, IncludeBatch)
{
    LR_model_formula <- switch(IncludeBatch,
            Yes = {
                ## mixed logistic model framework
                if(numberofsexes==2){
                    if (!sexIncluded  && !dimorphismIncluded){
                        as.formula(paste(depVariable, "~", paste("Genotype", "(1|Batch)", sep="+")))
                    } else if((sexIncluded && dimorphismIncluded)|(!sexIncluded && dimorphismIncluded)){
                        as.formula(paste(depVariable, "~",paste("Genotype", "Sex",  "Genotype*Sex", "(1|Batch)", sep= "+")))
                    } else if(sexIncluded  && !dimorphismIncluded){
                        as.formula(paste(depVariable, "~",paste("Genotype","Sex","(1|Batch)", sep= "+")))
                    }
                    
                }else{
                    as.formula(paste(depVariable, "~", paste("Genotype", "(1|Batch)", sep="+")))
                }
            },
            No = {
                ## standard logistic model framework
                if(numberofsexes==2){
                    if (!sexIncluded  && !dimorphismIncluded){
                        as.formula(paste(depVariable, "~", "Genotype", sep=""))
                    } else if((sexIncluded && dimorphismIncluded)|(!sexIncluded && dimorphismIncluded)){
                        as.formula(paste(depVariable, "~",paste("Genotype", "Sex",  "Genotype*Sex", sep= "+")))
                    } else if(sexIncluded  && !dimorphismIncluded){
                        as.formula(paste(depVariable, "~",paste("Genotype","Sex", sep= "+")))
                    }
                }else{
                    as.formula(paste(depVariable, "~","Genotype", sep=" "))
                }
            }
            )
    return(LR_model_formula)
}
##------------------------------------------------------------------------------
## Works with PhenTestResult object created by testDataset function.
## Builds null and uses final model to assess genotype effect
## final model results are then captured and stored.  
finalLRModel <- function(phenTestResult, outputMessages=TRUE)
{
    stop_message <- ""    
    ## START Check PhenTestResult object
    if(is(phenTestResult,"PhenTestResult")) {
        result <- phenTestResult
        x <- analysedDataset(result)
        linearRegressionOutput <- analysisResults(result)
        depVariable <- getVariable(result)
        numberofsexes <- linearRegressionOutput$numberSexes
        keep_sex <- linearRegressionOutput$model.effect.sex
        keep_interaction <- linearRegressionOutput$model.effect.interaction
        keep_batch <- linearRegressionOutput$model.effect.batch        
        ## Stop function if there are no datasets to work with
        if(is.null(x))
            stop_message <- "Error:\nPlease create a PhenList object first and run function 'testDataset'.\n"      
        ## Stop function if there are no enough input parameters
        if ( is.null(depVariable) || is.null(keep_batch)
                || is.null(keep_sex) || is.null(keep_interaction))
            stop_message <- "Error:\nPlease run function 'testDataset' first.\n"
    }
    else{
        stop_message <- "Error:\nPlease create a PhenTestResult object first.\n"
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
    ## END Checks and stop messages    
    result <- tryCatch( {   
                ## Test: genotype groups association with dependent variable
                ## Null Hypothesis: genotypes are not associated with dependent variable
                ## Alternative Hypothesis: genotypes are associated with dependent variable build genotype model
                # Formula is created for logistf method that doesn't support random effects 
                # so batch is excluded even if significant 
                if(numberofsexes==2){
                    if (!keep_sex  && !keep_interaction){
                        model_genotype.formula <- as.formula(paste(depVariable, "~", "Genotype", sep=""))
                    } else if((keep_sex && keep_interaction)|(!keep_sex && keep_interaction)){
                        model_genotype.formula <- as.formula(paste(depVariable, "~",paste("Sex",  "Genotype:Sex", sep= "+")))
                    } else if(keep_sex  && !keep_interaction){
                        model_genotype.formula <- as.formula(paste(depVariable, "~",paste("Genotype","Sex", sep= "+")))
                    }
                }else{
                    model_genotype.formula <- as.formula(paste(depVariable, "~","Genotype", sep=" "))
                }
               
                genotypeModel <- do.call("logistf", args = list(model_genotype.formula, data = x, na.action="na.omit" ))                        
                
                # Check if it is intercept only model    
                if(!((keep_sex && keep_interaction)|(!keep_sex && keep_interaction)|
                                (keep_sex  && !keep_interaction))){                                    
                    genotypeTest_p.value=logistftest(genotypeModel)$prob   
                    # !!!
                    # Alternate strategy needed to test when have intercept only null model as you cannot 
                    # specify an intercept only model with logistf       
                    linearRegressionOutput$model.null <- NA  
                    linearRegressionOutput$model.formula.null <- NA 
                }else{                    
                    # Build null model 
                    if(numberofsexes==2){
                        if (!keep_sex  && !keep_interaction){
                            model_null.formula <- as.formula(paste(depVariable, "~", 1, sep=""))
                        } else if((keep_sex && keep_interaction)|(!keep_sex && keep_interaction)|
                                (keep_sex  && !keep_interaction)){
                            model_null.formula <- as.formula(paste(depVariable, "~",paste( "Sex", sep= "+")))
                        }
                    }else{
                        model_null.formula <- as.formula(paste(depVariable, "~", "1", sep=" "))
                    }
                    nullModel <- do.call("logistf", args = list(model_null.formula, data = x, na.action="na.omit")) 
                    # Compare with genotype model
                    # genotypeTest_p.value <- pchisq((deviance(genotypeModel)-deviance(nullModel)), 1, lower=FALSE)  
                    # Problem with the estimate of df
                    genotypeTest_p.value = anova(genotypeModel, nullModel)$pval  
                    linearRegressionOutput$model.formula.null <- model_null.formula 
                    linearRegressionOutput$model.null <-nullModel                  
                }
                ## Store the results
                linearRegressionOutput$model.output <- genotypeModel
                linearRegressionOutput$model.formula.genotype <- model_genotype.formula
                linearRegressionOutput$model.output.genotype.nulltest.pVal <- genotypeTest_p.value                
                ## Create modeloutput and choose output depending on model
                linearRegressionOutput$model.output.summary  <- parserOutputSummaryLR(linearRegressionOutput)  
                result@analysisResults <- linearRegressionOutput 
                ## Assign LR quality of fit
                result@analysisResults$model.output.quality <- testFinalLRModel(result)
                result
            },            
            # End of tryCatch statement - if fails try to suggest something useful for the user
            error=function(error_mes) {
                result <- NULL
                stop_message <- paste("Error:\nCan't fit the model ", format(result$model_genotype.formula),
                        ". Try FE Framework.\n",sep="")
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
    )        
    return(result)
}
##------------------------------------------------------------------------------
## Creates model output summary and returns it in readable vector format
parserOutputSummaryLR<-function(linearRegressionOutput)
{
    # Estimates of coefficients
    coefficients <- format(linearRegressionOutput$model.output$coefficients,scientific=FALSE)
    # Standard errors
    error_estimates=format(sqrt(diag(vcov(linearRegressionOutput$model.output))),scientific=FALSE)
    # p-values
    probs <- format(linearRegressionOutput$model.output$prob,scientific=FALSE)        
    # Set all values to NA initially prior to selecting those relevant to the model 
    genotype_estimate <- NA
    genotype_estimate_SE <- NA
    genotype_p_value <- NA    
    sex_estimate <- NA
    sex_estimate_SE <- NA
    sex_p_value <- NA    
    intercept_estimate <- NA
    intercept_estimate_SE <- NA
    weight_estimate <- NA
    weight_estimate_SE <- NA
    weight_p_value <- NA    
    sex_FvKO_estimate <- NA
    sex_FvKO_SE <- NA
    sex_FvKO_p_value <- NA
    sex_MvKO_estimate <- NA
    sex_MvKO_SE <- NA
    sex_MvKO_p_value <- NA    
    # Pull out the information that is used to drive decision tree in where the information resides
    keep_sex <- linearRegressionOutput$model.effect.sex
    keep_interaction <- linearRegressionOutput$model.effect.interaction
    keep_batch <- linearRegressionOutput$model.effect.batch
    no_of_sexes <- linearRegressionOutput$numberSexes
    # Decision tree to pull the information depending on the final model fitted
    # note position is not dependent on whether a mixed or standard logisitic model is fitted.
    intercept_estimate <- coefficients[1]
    intercept_estimate_SE <- error_estimates[1]    
    if(no_of_sexes==1){        
        # Capturing SE based on advice in 
        # http://stats.stackexchange.com/questions/17571/
        # how-to-store-the-standard-errors-with-the-lm-function-in-r    
        genotype_estimate <- coefficients[2]
        genotype_estimate_SE <- error_estimates[2]
        genotype_p_value <- probs[2]        
    }else if(keep_interaction==TRUE){
        sex_estimate <- coefficients[2]
        sex_estimate_SE <- error_estimates[2]
        sex_p_value <- probs[2]    
        sex_FvKO_estimate <- coefficients[3]
        sex_FvKO_SE <- error_estimates[3]
        sex_FvKO_p_value <- probs[3]
        sex_MvKO_estimate <- coefficients[4]
        sex_MvKO_SE <- error_estimates[4]
        sex_MvKO_p_value <- probs[4]
    }else if(keep_interaction!=TRUE && keep_sex!=TRUE){
        genotype_estimate <- coefficients[2]
        genotype_estimate_SE <- error_estimates[2]
        genotype_p_value <- probs[2]
        sex_estimate <- coefficients[3]
        sex_estimate_SE <- error_estimates[3]
        sex_p_value <- probs[3]        
    }else{
        sex_estimate <- coefficients[2]
        sex_estimate_SE <- error_estimates[2]
        sex_p_value <- probs[2]
        genotype_estimate <- coefficients[3]
        genotype_estimate_SE <- error_estimates[3]
        genotype_p_value <- probs[3]    
    }    
    output <- c(genotype_estimate, genotype_estimate_SE, genotype_p_value, sex_estimate, sex_estimate_SE,  sex_p_value, 
            "NA", "NA", "NA", 
            intercept_estimate, intercept_estimate_SE, 
            sex_FvKO_estimate, sex_FvKO_SE, sex_FvKO_p_value,  
            sex_MvKO_estimate, sex_MvKO_SE, sex_MvKO_p_value)    
    names(output) <- c("genotype_estimate", "genotype_estimate_SE", 
            "genotype_p_value", 
            "sex_estimate", "sex_estimate_SE", "sex_p_value", 
            "weight_estimate", "weight_estimate_SE", "weight_p_value", 
            "intercept_estimate", "intercept_estimate_SE", 
            "sex_FvKO_estimate", "sex_FvKO_SE", "sex_FvKO_p_value", 
            "sex_MvKO_estimate", "sex_MvKO_SE", "sex_MvKO_p_value")
    return(output) 
}
##------------------------------------------------------------------------------