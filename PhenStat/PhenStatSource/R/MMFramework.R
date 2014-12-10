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
## MMFramework.R contains startModel, finalModel, modelFormula
## and parserOutputSummary functions
##------------------------------------------------------------------------------
## startModel creates start model and modifies it after testing of different
## hypothesis (the model effects).
## The model effects are:
## - batch effect (random effect significance),
## - variance effect (TRUE if residual variances for genotype groups are
## homogeneous and FALSE if they are heterogeneous),

## Fixed effects:
## - interaction effect (genotype by sex interaction significance),
## - sex effect (sex significance),
## - weigth effect (weigth significance).

## If user would like to assign other TRUE/FALSE values to the effects of
## the model then he or she has to define keepList argument
## which is vector of TRUE/FALSE values.

## If user has defined fixed effects (keepList argument) then function
## prints out calculated and user defined effects
## (only when outputMessages argument is set to TRUE),
## checks user defined effects for consistency
## (for instance, if there are no "Weight" column in the dataset then weight
## effect can't be assigned to TRUE, etc.)
## and modifies start model according to user defined effects.

## As the result PhenTestResult object that contains calculated or user
## defined model effects and MM start model is created.
startModel <- function(phenList, depVariable, equation="withWeight",
        outputMessages=TRUE, pThreshold=0.05, keepList=NULL)
{
    x <- phenList$dataset
    
    if (!is.null(keepList)){
        ## User's values for effects
        user_keep_weight <- keepList[3]
        user_keep_sex <- keepList[4]
        user_keep_interaction <- keepList[5]
        user_keep_batch <- keepList[1]
        user_keep_equalvar <- keepList[2]
        
        if (!('Weight' %in% colnames(x)) && user_keep_weight){
            if (outputMessages)
            message("Warning:\nWeight column is missed in the dataset. 'keepWeight' is set to FALSE.")
            
            user_keep_weight <- FALSE
        }
        
        if (!('Batch' %in% colnames(x)) && user_keep_batch){
            if (outputMessages)
            message("Warning:\nBatch column is missed in the dataset. 'keepBatch' is set to FALSE.")
            
            user_keep_batch <- FALSE
        }
        
    }
    
    numberofsexes <- length(levels(x$Sex))
    # Averages for percentage changes - is the ratio of the genotype effect for a sex relative to 
    # the wildtype signal for that variable for that sex - calculation        
    WT <- subset(x,x$Genotype==phenList$refGenotype)
    mean_all <- mean(WT[,c(depVariable)],na.rm=TRUE)  
    mean_list <- c(mean_all)  
    if (numberofsexes==2){  
        WT_f <- subset(WT,WT$Sex=="Female")
        WT_m <- subset(WT,WT$Sex=="Male")
        mean_f <- mean(WT_f[,c(depVariable)],na.rm=TRUE)
        mean_m <- mean(WT_m[,c(depVariable)],na.rm=TRUE)
        mean_list <- c(mean_all,mean_f,mean_m)  
    }
    # end of percentage change calculations    
    
    
    
    ## Start model formula: homogenous residual variance,
    ## genotype and sex interaction included
    model.formula  <- modelFormula(equation,numberofsexes, depVariable)
    
    #START OF tryCatch    
    finalResult <- tryCatch({
                
                if ('Batch' %in% colnames(x)){
                    ## GLS fit of model formula (no random effects)
                    ## Model 1A (model_withoutbatch)
                    model_withoutbatch <- do.call("gls",
                            args=list(model.formula, x, na.action="na.omit"))
                    ## MM fit of model formula (with random effects)
                    ## Model 1 (model_MM)
                    model_MM <- do.call("lme", args = list(model.formula,
                                            random=~1|Batch, data = x, na.action="na.omit", method="REML"))
                    ## Test: the random effects associated with batch intercepts can be
                    ## ommited from model
                    ## Hypothesis 1
                    ## Null Hypothesis: variance of batch = 0
                    ## Alternative Hypothesis: variance of batch > 0
                    ## For the division by 2 explanations see p.80 of "Linear Mixed Models"

                    p.value.batch <- (anova(model_MM, model_withoutbatch)$"p-value"[2])/2
                    ## The result of the test for Hypothesis 1 will help to select
                    ## the structure for random effects
                    keep_batch <- p.value.batch<pThreshold
                                        
                    ## MM fit of model formula with heterogeneous residual variances for
                    ## genotype groups
                    ## Model 1 assumes homogeneous residual variances
                    ## Model 2 with heterogeneous residual variances
                    model_hetvariance <- do.call("lme", args=list(model.formula,
                                            random=~1|Batch, x, weights=varIdent(form=~1|Genotype),
                                            na.action="na.omit", method="REML"))
                    

                    ## Test: the variance of the residuals is the same (homogeneous)
                    ## for all genotype groups
                    ## Hypothesis 2
                    ## Null Hypothesis: all residual variances are equal
                    ## Alternative Hypothesis: the residue variance is not equal
                    p.value.variance <- (anova(model_MM, model_hetvariance)$"p-value"[2])
                    ## The result of the test for Hypothesis 2 will help to select a
                    ## covariance structure for the residuals
                    keep_equalvar <- p.value.variance>pThreshold
                 
                }
                else {
                    ## No Batch effects
                    keep_batch <- FALSE
                    
                    ## Model 1A (model_withoutbatch)
                    model_MM <- do.call("gls",
                            args=list(model.formula, x, na.action="na.omit"))
                    
                    model_withoutbatch <- model_MM
                    
                    ## MM fit of model formula with heterogeneous residual variances for
                    ## genotype groups
                    ## Model 1 assumes homogeneous residual variances
                    ## Model 2 with heterogeneous residual variances
                    model_hetvariance <- do.call("gls", args=list(model.formula, x,
                                            weights=varIdent(form=~1|Genotype), na.action="na.omit"))

                    ## Test: the variance of the residuals is the same (homogeneous)
                    ## for all genotype groups
                    ## Hypothesis 2
                    ## Null Hypothesis: all residual variances are equal
                    ## Alternative Hypothesis: the residue variance is not equal
                    p.value.variance <- (anova(model_MM, model_hetvariance)$"p-value"[2])
                    ## The result of the test for Hypothesis 2 will help to select a
                    ## covariance structure for the residuals
                    keep_equalvar <- p.value.variance>pThreshold
                    
                    
                }
                
                
                ## Model fit is selected according to test results
                if(keep_batch && keep_equalvar){
                    ## Model 1
                    model <- model_MM
                }else if(keep_batch && !keep_equalvar){
                    ## Model 2
                    model <- model_hetvariance
                }else if(!keep_batch && keep_equalvar){
                    ## Model 1A
                    model= model_withoutbatch
                }else if(!keep_batch && !keep_equalvar){
                    ## Modify model 1A to heterogeneous residual variances
                    model <- do.call("gls", args=list(model.formula,
                                    weights=varIdent(form=~1|Genotype), x, na.action="na.omit"))
                }
                
                
                ## Tests for significance of fixed effects using TypeI F-test from anova
                ## functionality by using selected model
                anova_results <- anova(model, type="marginal")$"p-value" < pThreshold
                
                if(numberofsexes==2){
                    ## Result of the test for sex significance (fixed effect 1.)
                    keep_sex <- anova_results[3]
                    ## Eq.2
                    if (equation=="withWeight"){
                        ## Result of the test for weight significance  (fixed effect 3.)
                        
                        keep_weight <- anova_results[4]
                        ## Result of the test for genotype by sex interaction
                        ## significance (fixed effect 2.)
                        keep_interaction <- anova_results[5]
                        
                        ## Technical results needed for the output
                        ## Interaction test results are kept for the output
                        interactionTest <- anova(model, type="marginal")$"p-value"[5]
                        
                    }
                    ## Eq.1
                    else{
                        ## Result of the test for weight significance  (fixed effect 3.)
                        ## It's FALSE since here the equation 1 is used - without weight
                        ## effect
                        keep_weight <- FALSE
                        ## Result of the test for genotype by sex interaction
                        ## significance (fixed effect 2.)
                        keep_interaction <- anova_results[4]
                        ## Interaction test results are kept for the output
                        interactionTest <- anova(model, type="marginal")$"p-value"[4]
                    }
                }
                else {
                    keep_sex <- FALSE
                    keep_interaction <- FALSE
                    interactionTest <- NA
                    if (equation=="withWeight")
                    keep_weight <- anova_results[3]
                    else
                    keep_weight <- FALSE
                }
                
                if (!keep_weight && equation=="withWeight") {
                    equation="withoutWeight"
                    if (outputMessages)
                    message(paste("Since weight effect is not significant the equation ",
                                    "'withoutWeight' should be used instead.",sep=""))
                }
                
                if (outputMessages)
                message(paste("Information:\nEquation: '",equation,"'.\n",sep=""))
                
                if (outputMessages)
                message(paste("Information:\nCalculated values for model effects are: ",
                                "keepBatch=",keep_batch,
                                ", keepEqualVariance=",keep_equalvar,
                                ", keepWeight=",keep_weight,
                                ", keepSex=",keep_sex,
                                ", keepInteraction=",keep_interaction,".\n",sep=""))
                
                ## Results for user defined model effects values
                if (!is.null(keepList)){
                    if (outputMessages)
                    message(paste("Information:\nUser's values for model effects are: ",
                                    "keepBatch=",user_keep_batch,
                                    ", keepEqualVariance=",user_keep_equalvar,
                                    ", keepWeight=",user_keep_weight,
                                    ", keepSex=",user_keep_sex,
                                    ", keepInteraction=",user_keep_interaction,".\n",sep=""))
                    ## Model fit is selected according to user defined model effects
                    if(user_keep_batch && user_keep_equalvar){
                        ## Model 1
                        model <- model_MM
                    }else if(user_keep_batch && !user_keep_equalvar){
                        ## Model 2
                        model <- model_hetvariance
                    }else if(!user_keep_batch && user_keep_equalvar){
                        ## Model 1A
                        model <- model_withoutbatch
                    }else if(!user_keep_batch && !user_keep_equalvar){
                        ## Modify model 1A to heterogeneous residual variances
                        model <- do.call("gls", args=list(model.formula,
                                        weights=varIdent(form=~1|Genotype), x, na.action="na.omit"))
                    }
                    
                    if(numberofsexes==2){
                        if (equation=="withWeight"){
                            interactionTest <- anova(model, type="marginal")$"p-value"[5]
                            
                        }
                        else{
                            interactionTest <- anova(model, type="marginal")$"p-value"[4]
                            
                        }
                    }
                    else {
                        interactionTest <- NA
                    }
                    
                    compList <- (keepList==c(keep_batch,keep_equalvar,keep_weight,
                                    keep_sex,keep_interaction))
                    
                    if (length(compList[compList==FALSE])>0 && outputMessages)
                    message("Warning:\nCalculated values differ from user defined values for model effects.\n")
                    
                    keep_weight <- user_keep_weight
                    keep_sex <- user_keep_sex
                    keep_interaction <- user_keep_interaction
                    keep_batch <- user_keep_batch
                    keep_equalvar <- user_keep_equalvar
                    
                }
                
                
                finalResult <- new("PhenTestResult",list(
                                model.dataset=x,
                                model.output=model,
                                depVariable=depVariable,
                                refGenotype=phenList$refGenotype,
                                equation=equation,
                                method="MM",
                                model.effect.batch=keep_batch,
                                model.effect.variance=keep_equalvar,
                                model.effect.interaction=keep_interaction,
                                model.output.interaction=interactionTest,
                                model.effect.sex=keep_sex,
                                model.effect.weight=keep_weight,
                                numberSexes=numberofsexes,
                                pThreshold=pThreshold,
                                model.formula.genotype=model.formula,
                                model.output.averageRefGenotype=mean_list))
            },
            #END OF tryCatch    
            error=function(error_mes) {
                finalResult <- NULL
                if (equation=="withWeight") 
                stop_message <- paste("Error:\nCan't fit the model ",
                        format(model.formula),". Try MM with equation 'withoutWeight'. ",
                        "Another option is jitter.\n",sep="")
                else
                stop_message <- paste("Error:\nCan't fit the model ",
                        format(model.formula),". Try to add jitter or RR plus method.\n",sep="")
                
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
## Creates formula for the start model based on equation and number of sexes
## in the data
modelFormula <- function(equation, numberofsexes, depVariable)
{
    
    model.formula <- switch(equation,
            ## Eq.2
            withWeight = {
                ## Fixed effects: 1) Genotype 2) Sex 3) Genotype by Sex
                ## interaction 4) Weight
                if(numberofsexes==2){
                    as.formula(paste(depVariable, "~", paste("Genotype",
                                            "Sex", "Genotype*Sex","Weight", sep= "+")))
                }else{
                    as.formula(paste(depVariable, "~", paste("Genotype",
                                            "Weight",  sep= "+")))
                }
            },
            ## Eq.1
            withoutWeight = {
                ## Fixed effects: 1) Genotype 2) Sex 3) Genotype by Sex
                ## interaction
                if(numberofsexes==2){
                    as.formula(paste(depVariable, "~",
                                    paste("Genotype", "Sex", "Genotype*Sex", sep= "+")))
                }else{
                    as.formula(paste(depVariable, "~",
                                    paste("Genotype",  sep= "+")))
                }
            }
            )
    return(model.formula)
}

##------------------------------------------------------------------------------
## Works with PhenTestResult object created by testDataset function.
## Builds final model based on the significance of different model effects,
## depVariable and equation stored in phenTestResult object (see testDataset.R).
finalModel <- function(phenTestResult, outputMessages=TRUE)
{
    ## Checks and stop messages
    stop_message <- ""
    
    ## Check PhenTestResult object
    if(is(phenTestResult,"PhenTestResult")) {
        result <- phenTestResult
        x <- result$model.dataset
        depVariable <- result$depVariable
        equation <- result$equation
        keep_weight <- result$model.effect.weight
        keep_sex <- result$model.effect.sex
        keep_interaction <- result$model.effect.interaction
        keep_batch <- result$model.effect.batch
        keep_equalvar <- result$model.effect.variance
        
        ## Stop function if there are no datasets to work with
        if(is.null(x))
        stop_message <- "Error:\nPlease create a PhenList object first and run function 'testDataset'.\n"
        
        ## Stop function if there are no enough input parameters
        if (is.null(equation) || is.null(depVariable) || is.null(keep_batch)
                || is.null(keep_equalvar)
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
    
    numberofsexes <- result$numberSexes
    
    
    ## Build final null model
    ## Goal:  to test fixed effects of the model and based on the output
    ## build the final null model formula for later
    ## testing - as a null model it automatically excludes genotype and
    ## interaction term.
    ## The anova function tests the fixed effects associated by treatment
    ## with a null hypothesis that the regression
    ## coefficients are equal to zero  and an alternative hypothesis that the
    ## regression coefficient are not equal to zero.
    ## If the p-values of these tests are less than 0.05 we reject the null
    ## and accept the alternative that the are
    ## significant components of the model and should be included.
    ## If no terms are significant a model can be build with just an
    ## intercept element this is specified as
    ## "model.formula <- as.formula(paste(depVariable, "~", "1"))"
    
    ## Null model: genotype is not significant
    model_null.formula <- switch(equation,
            withWeight = {
                ## Eq.2
                if(numberofsexes==2){
                    if(!keep_sex){
                        as.formula(paste(depVariable, "~", "Weight"))
                        
                    }else{
                        as.formula(paste(depVariable, "~", 
                                        paste("Sex", "Weight", sep= "+")))
                    }
                }else{
                    as.formula(paste(depVariable, "~", "Weight"))
                }
            },
            withoutWeight = {
                ## Eq.1
                if(numberofsexes==2){
                    if(!keep_sex && !keep_interaction){
                        as.formula(paste(depVariable, "~", "1"))
                    }else{
                        as.formula(paste(depVariable, "~", "Sex"))
                    }
                }else{
                    as.formula(paste(depVariable, "~", "1"))
                }
            }
            )
    
    ## Alternative model: genotype is significant
    model_genotype.formula <- switch(equation,
            withWeight = {
                ## Eq.2
                if(numberofsexes==2){
                    if ((keep_sex && keep_weight && keep_interaction)|
                            (!keep_sex && keep_weight && keep_interaction)){
                        as.formula(paste(depVariable, "~",
                                        paste("Sex", "Genotype:Sex", "Weight", sep= "+")))
                        
                    } else if(keep_sex && keep_weight && !keep_interaction){
                        as.formula(paste(depVariable, "~",
                                        paste("Genotype", "Sex", "Weight", sep= "+")))
                        
                    } else if(!keep_sex && keep_weight && !keep_interaction){
                        as.formula(paste(depVariable, "~",
                                        paste("Genotype", "Weight", sep= "+")))
                    }
                    
                }else{
                    as.formula(paste(depVariable, "~",
                                    paste("Genotype", "Weight", sep="+")))
                }
            },
            withoutWeight = {
                ## Eq.1
                if(numberofsexes==2){
                    if (!keep_sex  && !keep_interaction){
                        as.formula(paste(depVariable, "~", "Genotype"))
                    } else if((keep_sex && keep_interaction)|
                            (!keep_sex && keep_interaction)){
                        as.formula(paste(depVariable, "~",
                                        paste("Sex",  "Genotype:Sex", sep= "+")))
                    } else if(keep_sex && !keep_interaction){
                        as.formula(paste(depVariable, "~",
                                        paste("Genotype","Sex", sep= "+")))
                    }
                }else{
                    as.formula(paste(depVariable, "~", paste("Genotype")))
                }
            }
            )
    
    finalResult <- tryCatch( {   
                ## Test: genotype groups association with dependent variable
                ## Null Hypothesis: genotypes are not associated with dependent variable
                ## Alternative Hypothesis: genotypes are associated with dependent
                ## variable
                if(keep_batch && keep_equalvar){

                    model_genotype <-  do.call("lme", args = list(model_genotype.formula,
                                            random=~1|Batch, x, na.action="na.omit", method="ML"))
                    
                    model_null <- do.call("lme", args=list(model_null.formula, x,
                                    random=~1|Batch, na.action="na.omit",  method="ML"))
                    p.value <- (anova(model_genotype, model_null)$"p-value"[2])
                   
                }
                if(keep_batch && !keep_equalvar){
                    model_genotype <- do.call("lme", args = list(model_genotype.formula,
                                    random=~1|Batch, x,weights=varIdent(form=~1|Genotype),
                                    na.action="na.omit", method="ML"))
                    
                    model_null <- do.call("lme", args=list(model_null.formula, x,
                                    random=~1|Batch,weights=varIdent(form=~1|Genotype),
                                    na.action="na.omit",  method="ML"))
                    
                    p.value <- (anova(model_genotype, model_null)$"p-value"[2])
                }else if(!keep_batch && !keep_equalvar){
                    model_genotype <- do.call("gls", args = list(model_genotype.formula,
                                    x,weights=varIdent(form=~1|Genotype),method="ML", na.action="na.omit"))
                    
                    model_null <- do.call("gls", args=list(model_null.formula,
                                    x,weights=varIdent(form=~1|Genotype),method="ML", na.action="na.omit"))
                    
                    p.value <- (anova(model_genotype, model_null)$"p-value"[2])
                }else if(!keep_batch && keep_equalvar){
                    model_genotype <- do.call("gls", args = list(model_genotype.formula,
                                    x, method="ML", na.action="na.omit"))
                    
                    model_null <- do.call("gls", args=list(model_null.formula, x,
                                    method="ML", na.action="na.omit"))
                    
                    p.value <- (anova(model_genotype, model_null)$"p-value"[2])
                }
                
                ## Final model version with na.exclude and REML method
                if(keep_batch && keep_equalvar){
                    ## Model 1
                    model_genotype <- do.call("lme", args = list(model_genotype.formula,
                                    random=~1|Batch, x, na.action="na.exclude", method="REML"))
                }else if(keep_batch && !keep_equalvar){
                    ## Model 2
                    model_genotype <- do.call("lme", args = list(model_genotype.formula,
                                    random=~1|Batch, x,weights=varIdent(form=~1|Genotype),
                                    na.action="na.exclude", method="REML"))
                }else if(!keep_batch && !keep_equalvar){
                    ## Model 2A
                    model_genotype <- do.call("gls", args = list(model_genotype.formula,
                                    x,weights=varIdent(form=~1|Genotype), na.action="na.exclude"))
                }else if(!keep_batch && keep_equalvar){
                    ## Model 1A
                    model_genotype <- do.call("gls", args = list(model_genotype.formula,
                                    x, na.action="na.exclude"))
                }
                
                
                ## Store the results
                result$model.output <- model_genotype
                result$model.null <- model_null
                result$model.output.genotype.nulltest.pVal <- p.value
                result$model.formula.null <- model_null.formula
                result$model.formula.genotype <- model_genotype.formula
                result$model.effect.variance <- keep_equalvar
                
                ## Assign MM quality of fit
                result$model.output.quality <- testFinalModel(result)
                
                ## Parse modeloutput and choose output depending on model
                result$model.output.summary <- parserOutputSummary(result)

                # Percentage changes - is the ratio of the genotype effect for a sex relative to 
                # the wildtype signal for that variable for that sex - calculation     

                if(result$numberSexes==2){
                   # without weight
                   if (is.na(result$model.output.summary['weight_estimate'])){ 
                        if (!is.na(result$model.output.summary['sex_estimate']) &&
                            !is.na(result$model.output.summary['sex_FvKO_estimate']))
                        {
                            denominator_f <- result$model.output.summary['intercept_estimate']
                            denominator_m <- result$model.output.summary['intercept_estimate']+
                            result$model.output.summary['sex_estimate']
                            ratio_f <- result$model.output.summary['sex_FvKO_estimate']/denominator_f                       
                            ratio_m <- result$model.output.summary['sex_MvKO_estimate']/denominator_m
                        }
                        else if (!is.na(result$model.output.summary['sex_FvKO_estimate']))
                        {
                                denominator <- result$model.output.summary['intercept_estimate']
                                ratio_f <- result$model.output.summary['sex_FvKO_estimate']/denominator                            
                                ratio_m <- result$model.output.summary['sex_MvKO_estimate']/denominator            
                        }
                        else if (!is.na(result$model.output.summary['sex_estimate']))
                        {
                            denominator_f <- result$model.output.summary['intercept_estimate']
                            denominator_m <- result$model.output.summary['intercept_estimate']+
                            result$model.output.summary['sex_estimate']
                            ratio_f <- result$model.output.summary['genotype_estimate']/denominator_f                        
                            ratio_m <- result$model.output.summary['genotype_estimate']/denominator_m 
                        }
                        else
                        {
                            denominator <- result$model.output.summary['intercept_estimate']
                            ratio_f <- result$model.output.summary['genotype_estimate']/denominator                            
                            ratio_m <- ratio_f                      
                        }
                   }
                   # with weight
                   else{
                        mean_list <- result$model.output.averageRefGenotype
                        denominator_f <- mean_list[2]
                        denominator_m <- mean_list[3]
                        if (!is.na(result$model.output.summary['sex_estimate']) &&
                        !is.na(result$model.output.summary['sex_FvKO_estimate']))
                        {
                            ratio_f <- result$model.output.summary['sex_FvKO_estimate']/denominator_f                       
                            ratio_m <- result$model.output.summary['sex_MvKO_estimate']/denominator_m 
                        }
                        else if (!is.na(result$model.output.summary['sex_FvKO_estimate']))
                        {
                            ratio_f <- result$model.output.summary['sex_FvKO_estimate']/denominator_f                            
                            ratio_m <- result$model.output.summary['sex_MvKO_estimate']/denominator_m            
                        }
                        else 
                        {
                            ratio_f <- result$model.output.summary['genotype_estimate']/denominator_f                        
                            ratio_m <- result$model.output.summary['genotype_estimate']/denominator_m 
                        }
                   }

                   result$model.output.percentageChanges <- c(ratio_f*100,ratio_m*100)
                   names(result$model.output.percentageChanges) <- c('female*genotype ratio','male*genotype ratio')
                   finalResult <- result
                }
                else{
                    # without weight
                    if (is.na(result$model.output.summary['weight_estimate'])){ 
                        denominator <- result$model.output.summary['intercept_estimate']
                        ratio_f <- result$model.output.summary['genotype_estimate']/denominator                            
                    }
                    # with weight
                    else{
                        mean_list <- result$model.output.averageRefGenotype
                        denominator <- mean_list[1]
                        ratio_f <- result$model.output.summary['genotype_estimate']/denominator   
                    }

                    result$model.output.percentageChanges <- c(ratio_f*100)
                    names(result$model.output.percentageChanges) <- c('all*genotype ratio')
                    finalResult <- result
                }
                # end of percentage changes calculation


            },
            
            # End of tryCatch statement - if fails try to suggest smth useful for the user
            error=function(error_mes) {
                
                finalResult <- NULL
                if (equation=="withWeight") 
                stop_message <- paste("Error:\nCan't fit the model ",
                        format(model_genotype.formula),". Try MM with equation 'withoutWeight'. ",
                        "Another option is jitter\n",sep="")
                else
                stop_message <- paste("Error:\nCan't fit the model ",
                        format(model_genotype.formula),". Try to add jitter or RR plus method.\n",sep="")
                
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
## Parses model output summary and returns in readable vector format
parserOutputSummary<-function(phenTestResult)
{
    result <- phenTestResult
    modeloutput_summary <- summary(result$model.output)
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
    
    lengthoftable <- {
        table_length <- NA
        
        if (result$equation=="withWeight"){
            if(result$numberSexes==2){
                if((result$model.effect.sex
                                && result$model.effect.interaction)|
                        (!result$model.effect.sex && result$model.effect.interaction)){
                    table_length <- 5
                }else if(result$model.effect.sex &&
                        !result$model.effect.interaction){
                    table_length <- 4
                }else {
                    table_length <- 3
                }
            }else{
                table_length <- 3
            }
        }
        ## Eq.1
        else{
            if(result$numberSexes==2){
                if((result$model.effect.sex
                                && result$model.effect.interaction)|
                        (!result$model.effect.sex && result$model.effect.interaction)){
                    table_length <- 4
                }else if(!result$model.effect.sex &&
                        !result$model.effect.interaction){
                    table_length <- 2
                }else{
                    table_length <- 3
                }
            }else{
                table_length <- 2
            }
        }
        
        table_length
    }
    
    switch(result$equation,
            withoutWeight = {
                sex_index <- match(c("SexMale"),row.names(modeloutput_summary[["tTable"]]))
                sex_FvKO_index <- 3
                sex_MvKO_index <- 4
                if (is.na(sex_index)){
                    sex_index <- match(c("SexFemale"),row.names(modeloutput_summary[["tTable"]]))
                    sex_FvKO_index <- 4
                    sex_MvKO_index <- 3
                }
                if(result$model.effect.batch){
                    ## for mixed model
                    intercept_estimate = modeloutput_summary[["tTable"]][[1]]
                    intercept_estimate_SE = modeloutput_summary[["tTable"]][[(1+lengthoftable)]]
                    if((result$model.effect.sex && result$model.effect.interaction)
                            |( !result$model.effect.sex&& result$model.effect.interaction)){
                        sex_estimate=modeloutput_summary[["tTable"]][[sex_index]]
                        sex_estimate_SE=modeloutput_summary[["tTable"]][[(sex_index+lengthoftable)]]
                        sex_p_value= modeloutput_summary[["tTable"]][[(sex_index+4*lengthoftable)]]
                        sex_FvKO_estimate= modeloutput_summary[["tTable"]][[sex_FvKO_index]]
                        sex_FvKO_SE=modeloutput_summary[["tTable"]][[(sex_FvKO_index+lengthoftable)]]
                        sex_FvKO_p_value=modeloutput_summary[["tTable"]][[(sex_FvKO_index+4*lengthoftable)]]
                        sex_MvKO_estimate=modeloutput_summary[["tTable"]][[sex_MvKO_index]]
                        sex_MvKO_SE=modeloutput_summary[["tTable"]][[(sex_MvKO_index+lengthoftable)]]
                        sex_MvKO_p_value=modeloutput_summary[["tTable"]][[(sex_MvKO_index+4*lengthoftable)]]
                    } else if( !result$model.effect.sex && !result$model.effect.interaction){
                        genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                        genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                        genotype_p_value =  modeloutput_summary[["tTable"]][[(2+4*lengthoftable)]]
                        
                    }else{
                        genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                        genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                        genotype_p_value =  modeloutput_summary[["tTable"]][[(2+4*lengthoftable)]]
                        sex_estimate=modeloutput_summary[["tTable"]][[sex_index]]
                        sex_estimate_SE=modeloutput_summary[["tTable"]][[(sex_index+lengthoftable)]]
                        sex_p_value= modeloutput_summary[["tTable"]][[(sex_index+4*lengthoftable)]]
                    }
                    
                }else{
                    ## adaption for being a linear model rather than a mixed
                    ## model
                    intercept_estimate = modeloutput_summary[["tTable"]][[1]]
                    intercept_estimate_SE = modeloutput_summary[["tTable"]][[(1+lengthoftable)]]
                    
                    if((result$model.effect.sex && result$model.effect.interaction)
                            |( !result$model.effect.sex&& result$model.effect.interaction)){
                        
                        sex_estimate=modeloutput_summary[["tTable"]][[sex_index]]
                        sex_estimate_SE=modeloutput_summary[["tTable"]][[(sex_index+lengthoftable)]]
                        sex_p_value= modeloutput_summary[["tTable"]][[(sex_index+3*lengthoftable)]]
                        sex_FvKO_estimate= modeloutput_summary[["tTable"]][[sex_FvKO_index]]
                        sex_FvKO_SE=modeloutput_summary[["tTable"]][[(sex_FvKO_index+lengthoftable)]]
                        sex_FvKO_p_value=modeloutput_summary[["tTable"]][[(sex_FvKO_index+3*lengthoftable)]]
                        sex_MvKO_estimate=modeloutput_summary[["tTable"]][[sex_MvKO_index]]
                        sex_MvKO_SE=modeloutput_summary[["tTable"]][[(sex_MvKO_index+lengthoftable)]]
                        sex_MvKO_p_value=modeloutput_summary[["tTable"]][[(sex_MvKO_index+3*lengthoftable)]]
                        
                        
                    } else if( !result$model.effect.sex && !result$model.effect.interaction){
                        
                        genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                        genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                        genotype_p_value =  modeloutput_summary[["tTable"]][[(2+3*lengthoftable)]]
                        
                    }else{
                        
                        genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                        genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                        genotype_p_value =  modeloutput_summary[["tTable"]][[(2+3*lengthoftable)]]
                        sex_estimate=modeloutput_summary[["tTable"]][[sex_index]]
                        sex_estimate_SE=modeloutput_summary[["tTable"]][[(sex_index+lengthoftable)]]
                        sex_p_value= modeloutput_summary[["tTable"]][[(3+3*lengthoftable)]]
                    }
                    
                    
                }
            },
            withWeight = {
                sex_index <- match(c("SexMale"),row.names(modeloutput_summary[["tTable"]]))
                sex_FvKO_index <- 4
                sex_MvKO_index <- 5
                if (is.na(sex_index)){
                    sex_index <- match(c("SexFemale"),row.names(modeloutput_summary[["tTable"]]))
                    sex_FvKO_index <- 5
                    sex_MvKO_index <- 4
                }
                if(!result$model.effect.weight){
                    
                    ## If weight is not significant then the output is the
                    ## same as fitting model Eq1 and so no output is needed.
                    result$model.effect.batch=NA
                    
                }else{
                    
                    if(result$model.effect.batch){
                        
                        ## for mixed model
                        intercept_estimate = modeloutput_summary[["tTable"]][[1]]
                        intercept_estimate_SE = modeloutput_summary[["tTable"]][[(1+lengthoftable)]]
                        
                        if((result$model.effect.weight && result$model.effect.sex &&
                                        result$model.effect.interaction) |
                                (result$model.effect.weight &&
                                        !result$model.effect.sex &&
                                        result$model.effect.interaction)){
                            sex_estimate=modeloutput_summary[["tTable"]][[sex_index]]
                            sex_estimate_SE=modeloutput_summary[["tTable"]][[(sex_index+lengthoftable)]]
                            sex_p_value= modeloutput_summary[["tTable"]][[(sex_index+4*lengthoftable)]]
                            sex_FvKO_estimate= modeloutput_summary[["tTable"]][[sex_FvKO_index]]
                            sex_FvKO_SE=modeloutput_summary[["tTable"]][[(sex_FvKO_index+lengthoftable)]]
                            sex_FvKO_p_value=modeloutput_summary[["tTable"]][[(sex_FvKO_index+4*lengthoftable)]]
                            sex_MvKO_estimate=modeloutput_summary[["tTable"]][[sex_MvKO_index]]
                            sex_MvKO_SE=modeloutput_summary[["tTable"]][[(sex_MvKO_index+lengthoftable)]]
                            sex_MvKO_p_value=modeloutput_summary[["tTable"]][[(sex_MvKO_index+4*lengthoftable)]]
                            weight_estimate=modeloutput_summary[["tTable"]][[3]]
                            weight_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                            weight_p_value=modeloutput_summary[["tTable"]][[(3+4*lengthoftable)]]
                            
                        } else if (result$model.effect.weight &&
                                !result$model.effect.sex &&
                                !result$model.effect.interaction){
                            
                            genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput_summary[["tTable"]][[(2+4*lengthoftable)]]
                            weight_estimate=modeloutput_summary[["tTable"]][[3]]
                            weight_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                            weight_p_value=modeloutput_summary[["tTable"]][[(3+4*lengthoftable)]]
                            
                        }else if (result$model.effect.weight &&
                                result$model.effect.sex &&
                                !result$model.effect.interaction){
                            
                            genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput_summary[["tTable"]][[(2+4*lengthoftable)]]
                            sex_estimate=modeloutput_summary[["tTable"]][[3]]
                            sex_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                            sex_p_value= modeloutput_summary[["tTable"]][[(3+4*lengthoftable)]]
                            weight_estimate=modeloutput_summary[["tTable"]][[4]]
                            weight_estimate_SE=modeloutput_summary[["tTable"]][[(4+lengthoftable)]]
                            weight_p_value=modeloutput_summary[["tTable"]][[(4+4*lengthoftable)]]
                            
                        }
                        
                    }else{
                        ## adaption for being a linear model rather than a
                        ## mixed model
                        intercept_estimate = modeloutput_summary[["tTable"]][[1]]
                        intercept_estimate_SE = modeloutput_summary[["tTable"]][[(1+lengthoftable)]]
                        
                        if((result$model.effect.weight &&
                                        result$model.effect.sex &&
                                        result$model.effect.interaction )|
                                (result$model.effect.weight &&
                                        !result$model.effect.sex &&
                                        result$model.effect.interaction)){
                            sex_estimate=modeloutput_summary[["tTable"]][[sex_index]]
                            sex_estimate_SE=modeloutput_summary[["tTable"]][[(sex_index+lengthoftable)]]
                            sex_p_value= modeloutput_summary[["tTable"]][[(sex_index+3*lengthoftable)]]
                            weight_estimate=modeloutput_summary[["tTable"]][[3]]
                            weight_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                            weight_p_value=modeloutput_summary[["tTable"]][[(3+3*lengthoftable)]]
                            sex_FvKO_estimate= modeloutput_summary[["tTable"]][[sex_FvKO_index]]
                            sex_FvKO_SE=modeloutput_summary[["tTable"]][[(sex_FvKO_index+lengthoftable)]]
                            sex_FvKO_p_value=modeloutput_summary[["tTable"]][[(sex_FvKO_index+3*lengthoftable)]]
                            sex_MvKO_estimate=modeloutput_summary[["tTable"]][[sex_MvKO_index]]
                            sex_MvKO_SE=modeloutput_summary[["tTable"]][[(sex_MvKO_index+lengthoftable)]]
                            sex_MvKO_p_value=modeloutput_summary[["tTable"]][[(sex_MvKO_index+3*lengthoftable)]]
                            
                        } else if (result$model.effect.weight &&
                                result$model.effect.sex &&
                                !result$model.effect.interaction){
                            
                            genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput_summary[["tTable"]][[(2+3*lengthoftable)]]
                            sex_estimate=modeloutput_summary[["tTable"]][[sex_index]]
                            sex_estimate_SE=modeloutput_summary[["tTable"]][[(sex_index+lengthoftable)]]
                            sex_p_value= modeloutput_summary[["tTable"]][[(sex_index+3*lengthoftable)]]
                            weight_estimate=modeloutput_summary[["tTable"]][[4]]
                            weight_estimate_SE=modeloutput_summary[["tTable"]][[(4+lengthoftable)]]
                            weight_p_value=modeloutput_summary[["tTable"]][[(4+3*lengthoftable)]]
                            
                            
                        }else if (result$model.effect.weight &&
                                !result$model.effect.sex &&
                                !result$model.effect.interaction){
                            genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput_summary[["tTable"]][[(2+3*lengthoftable)]]
                            weight_estimate=modeloutput_summary[["tTable"]][[3]]
                            weight_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                            weight_p_value=modeloutput_summary[["tTable"]][[(3+3*lengthoftable)]]
                        }
                        
                    }
                    
                }
            }
            )
    output <- c(genotype_estimate, genotype_estimate_SE,  genotype_p_value,
            sex_estimate, sex_estimate_SE,  sex_p_value, 
            weight_estimate, weight_estimate_SE, weight_p_value, 
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