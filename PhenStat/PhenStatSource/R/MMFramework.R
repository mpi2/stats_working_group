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
## - interaction effect (genotype by gender interaction significance),
## - gender effect (gender significance),
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
        user_keep_gender <- keepList[4]
        user_keep_interaction <- keepList[5]
        user_keep_batch <- keepList[1]
        user_keep_equalvar <- keepList[2]
        
        if (!('Weight' %in% colnames(x)) && user_keep_weight){
            if (outputMessages)
            message("Warning:\nWeight column is missed in the dataset.
                    'keepWeight' is set to FALSE.")
            
            user_keep_weight <- FALSE
        }
        
        if (!('Batch' %in% colnames(x)) && user_keep_batch){
            if (outputMessages)
            message("Warning:\nBatch column is missed in the dataset.
                    'keepBatch' is set to FALSE.")
            
            user_keep_batch <- FALSE
        }
        
    }
    
    numberofgenders <- length(levels(x$Gender))
    
    
    ## Start model formula: homogenous residual variance,
    ## genotype and sex interaction included
    model.formula  <- modelFormula(equation,numberofgenders, depVariable)
    
    if ('Batch' %in% colnames(x)){
        ## GLS fit of model formula (no random effects)
        ## Model 1A (model_withoutbatch)
        model_withoutbatch <- do.call("gls",
                args=list(model.formula, x, na.action="na.omit"))
        ## MM fit of model formula (with random effects)
        ## Model 1 (model_MM)
        model_MM <-
        tryCatch(
                model_MM <- do.call("lme", args = list(model.formula,
                                random=~1|Batch, data = x, na.action="na.omit", method="REML")),
                error=function(error_mes) {
                    if (outputMessages)
                    message("Warning:\nMixed model with batch effect
                            as random effect is not fitting - false convergence.\n
                            Mixed model with no random effects is used instead.\n")
                    
                    model_MM <- NULL
                }
                )
        ## Test: the random effects associated with batch intercepts can be
        ## ommited from model
        ## Hypothesis 1
        ## Null Hypothesis: variance of batch = 0
        ## Alternative Hypothesis: variance of batch > 0
        ## For the division by 2 explanations see p.80 of "Linear Mixed Models"
        if (!is.null(model_MM)) {
            p.value.batch <- (anova(model_MM, model_withoutbatch)$p[2])/2
            ## The result of the test for Hypothesis 1 will help to select
            # the structure for random effects
            keep_batch <- p.value.batch<pThreshold
            
        }
        else {
            keep_batch <- FALSE
            if (!is.null(keepList)){
                if (outputMessages && user_defined_batch)
                message("Warning:\n'keepBatch' is set to FALSE otherwise
                        the model can't be fitted - false convergence.\n")
                
                user_defined_batch <- FALSE
            }
            model_MM <- model_withoutbatch
        }
        
        ## MM fit of model formula with heterogeneous residual variances for
        ## genotype groups
        ## Model 1 assumes homogeneous residual variances
        ## Model 2 with heterogeneous residual variances
        model_hetvariance <-
        tryCatch(
                model_hetvariance <- do.call("lme", args=list(model.formula,
                                random=~1|Batch, x, weights=varIdent(form=~1|Genotype),
                                na.action="na.omit", method="REML")),
                error=function(error_mes) {
                    if (outputMessages)
                    message("Warning:\nMixed model with heterogeneous
                            residual variances for genotype groups is not
                            fitting - false convergence.\nMixed model with
                            homogeneous residual variances is used instead.\n")
                    
                    model_hetvariance <- NULL
                }
                )
        
        if (!is.null(model_hetvariance)) {
            ## Test: the variance of the residuals is the same (homogeneous)
            ## for all genotype groups
            ## Hypothesis 2
            ## Null Hypothesis: all residual variances are equal
            ## Alternative Hypothesis: the residue variance is not equal
            p.value.variance <- (anova(model_MM, model_hetvariance)$p[2])
            ## The result of the test for Hypothesis 2 will help to select a
            ## covariance structure for the residuals
            keep_equalvar <- p.value.variance>pThreshold
        }
        else {
            keep_equalvar <- TRUE
            if (!is.null(keepList)){
                if (outputMessages && !user_keep_equalvar)
                message("Warning:\n'keepEqualVariance' is set to TRUE
                        otherwise the model can't be fitted - false convergence.\n")
                
                user_keep_equalvar <- TRUE
            }
        }
    }
    else {
        ## No Batch effects
        keep_batch <- FALSE
        
        ## Model 1A (model_withoutbatch)
        model_MM <- do.call("gls",
                args=list(model.formula, x, na.action="na.omit"))
        
        ## MM fit of model formula with heterogeneous residual variances for
        ## genotype groups
        ## Model 1 assumes homogeneous residual variances
        ## Model 2 with heterogeneous residual variances
        model_hetvariance <-
        tryCatch(
                model_hetvariance <- do.call("gls", args=list(model.formula, x,
                                weights=varIdent(form=~1|Genotype), na.action="na.omit")),
                error=function(error_mes) {
                    if (outputMessages)
                    message("Warning:\nMixed model with heterogeneous
                            residual variances for genotype groups is not
                            fitting - false convergence.\nMixed model with
                            homogeneous residual variances is used instead.\n")
                    
                    model_hetvariance <- NULL
                }
                )
        
        if (!is.null(model_hetvariance)) {
            ## Test: the variance of the residuals is the same (homogeneous)
            ## for all genotype groups
            ## Hypothesis 2
            ## Null Hypothesis: all residual variances are equal
            ## Alternative Hypothesis: the residue variance is not equal
            p.value.variance <- (anova(model_MM, model_hetvariance)$p[2])
            ## The result of the test for Hypothesis 2 will help to select a
            ## covariance structure for the residuals
            keep_equalvar <- p.value.variance>pThreshold
        }
        else {
            keep_equalvar <- TRUE
            if (!is.null(keepList)){
                if (outputMessages && !user_keep_equalvar)
                message("Warning:\n'keepEqualVariance' is set to TRUE
                        otherwise the model can't be fitted - false convergence.\n")
                
                user_keep_equalvar <- TRUE
            }
        }
        
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
    
    if(numberofgenders==2){
        ## Result of the test for gender significance (fixed effect 1.)
        keep_gender <- anova_results[3]
        ## Eq.2
        if (equation=="withWeight"){
            ## Result of the test for weight significance  (fixed effect 3.)
            
            keep_weight <- anova_results[4]
            ## Result of the test for genotype by gender interaction
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
            ## Result of the test for genotype by gender interaction
            ## significance (fixed effect 2.)
            keep_interaction <- anova_results[4]
            ## Interaction test results are kept for the output
            interactionTest <- anova(model, type="marginal")$"p-value"[4]
        }
    }
    else {
        keep_gender <- FALSE
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
        message("Since weight effect is not significant the equation
                Eq.1 'withoutWeight' should be used instead.")
    }
    
    if (outputMessages)
    message(paste("Information:\nCalculated values for model effects are:
                    keepBatch=",keep_batch,
                    ", keepEqualVariance=",keep_equalvar,
                    ", keepWeight=",keep_weight,
                    ", keepGender=",keep_gender,
                    ", keepInteraction=",keep_interaction,".\n",sep=""))
    
    ## Results for user defined model effects values
    if (!is.null(keepList)){
        if (outputMessages)
        message(paste("Information:\nUser's values for model effects are:
                        keepBatch=",user_keep_batch,
                        ", keepEqualVariance=",user_keep_equalvar,
                        ", keepWeight=",user_keep_weight,
                        ", keepGender=",user_keep_gender,
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
        
        if(numberofgenders==2){
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
                        keep_gender,keep_interaction))
        
        if (length(compList[compList==FALSE])>0 && outputMessages)
        message("Warning:\nCalculated values differ from user defined values
                for model effects.\n")
        
        keep_weight <- user_keep_weight
        keep_gender <- user_keep_gender
        keep_interaction <- user_keep_interaction
        keep_batch <- user_keep_batch
        keep_equalvar <- user_keep_equalvar
        
    }
    
    if (outputMessages)
    message(paste("Information:\nEquation: '",equation,"'.\n",sep=""))
    
    
    result <- new("PhenTestResult",list(
                    model.dataset=x,
                    model.output=model,
                    depVariable=depVariable,
                    equation=equation,
                    method="MM",
                    model.effect.batch=keep_batch,
                    model.effect.variance=keep_equalvar,
                    model.effect.interaction=keep_interaction,
                    model.output.interaction=interactionTest,
                    model.effect.gender=keep_gender,
                    model.effect.weight=keep_weight,
                    numberGenders=numberofgenders,
                    pThreshold=pThreshold,
                    model.formula.genotype=model.formula))
    return(result)
}

##------------------------------------------------------------------------------
## Creates formula for the start model based on equation and number of genders
## in the data
modelFormula <- function(equation, numberofgenders, depVariable)
{
    
    model.formula <- switch(equation,
            ## Eq.2
            withWeight = {
                ## Fixed effects: 1) Genotype 2) Gender 3) Genotype by Gender
                ## interaction 4) Weight
                if(numberofgenders==2){
                    as.formula(paste(depVariable, "~", paste("Genotype",
                                            "Gender", "Genotype*Gender","Weight", sep= "+")))
                }else{
                    as.formula(paste(depVariable, "~", paste("Genotype",
                                            "Weight",  sep= "+")))
                }
            },
            ## Eq.1
            withoutWeight = {
                ## Fixed effects: 1) Genotype 2) Gender 3) Genotype by Gender
                ## interaction
                if(numberofgenders==2){
                    as.formula(paste(depVariable, "~",
                                    paste("Genotype", "Gender", "Genotype*Gender", sep= "+")))
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
        keep_gender <- result$model.effect.gender
        keep_interaction <- result$model.effect.interaction
        keep_batch <- result$model.effect.batch
        keep_equalvar <- result$model.effect.variance
        
        ## Stop function if there are no datasets to work with
        if(is.null(x))
        stop_message <- "Error:\nPlease create a PhenList object first and
        run function 'testDataset'.\n"
        
        ## Stop function if there are no enough input parameters
        if (is.null(equation) || is.null(depVariable) || is.null(keep_batch)
                || is.null(keep_equalvar)
                || is.null(keep_gender) || is.null(keep_interaction))
        stop_message <- "Error:\nPlease run function 'testDataset' first.\n"
    }
    else{
        stop_message <- "Error:\nPlease create a PhenTestResult object first.\n"
    }
    
    
    if (nchar(stop_message)>1){
        if (outputMessages)
        message(stop_message)
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))
        stop()
    }
    
    ## END Checks and stop messages
    
    numberofgenders <- result$numberGenders
    
    
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
                if(numberofgenders==2){
                    if(!keep_gender){
                        as.formula(paste(depVariable, "~", "Weight"))
                        
                    }else{
                        as.formula(paste(depVariable, "~", 
                                        paste("Gender", "Weight", sep= "+")))
                    }
                }else{
                    as.formula(paste(depVariable, "~", "Weight"))
                }
            },
            withoutWeight = {
                ## Eq.1
                if(numberofgenders==2){
                    if(!keep_gender && !keep_interaction){
                        as.formula(paste(depVariable, "~", "1"))
                    }else{
                        as.formula(paste(depVariable, "~", "Gender"))
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
                if(numberofgenders==2){
                    if ((keep_gender && keep_weight && keep_interaction)|
                            (!keep_gender && keep_weight && keep_interaction)){
                        as.formula(paste(depVariable, "~",
                                        paste("Gender", "Genotype:Gender", "Weight", sep= "+")))
                        
                    } else if(keep_gender && keep_weight && !keep_interaction){
                        as.formula(paste(depVariable, "~",
                                        paste("Genotype", "Gender", "Weight", sep= "+")))
                        
                    } else if(!keep_gender && keep_weight && !keep_interaction){
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
                if(numberofgenders==2){
                    if (!keep_gender  && !keep_interaction){
                        as.formula(paste(depVariable, "~", "Genotype"))
                    } else if((keep_gender && keep_interaction)|
                            (!keep_gender && keep_interaction)){
                        as.formula(paste(depVariable, "~",
                                        paste("Gender",  "Genotype:Gender", sep= "+")))
                    } else if(keep_gender && !keep_interaction){
                        as.formula(paste(depVariable, "~",
                                        paste("Genotype","Gender", sep= "+")))
                    }
                }else{
                    as.formula(paste(depVariable, "~", paste("Genotype")))
                }
            }
            )
    
    
    ## Test: genotype groups association with dependent variable
    ## Null Hypothesis: genotypes are not associated with dependent variable
    ## Alternative Hypothesis: genotypes are associated with dependent
    ## variable
    if(keep_batch && keep_equalvar){
        model_genotype <- do.call("lme", args = list(model_genotype.formula,
                        random=~1|Batch, x, na.action="na.omit", method="ML"))
        
        model_null <- do.call("lme", args=list(model_null.formula, x,
                        random=~1|Batch, na.action="na.omit",  method="ML"))
        
        p.value <- (anova(model_genotype, model_null)$p[2])
    }else if(keep_batch && !keep_equalvar){
        model_genotype <- do.call("lme", args = list(model_genotype.formula,
                        random=~1|Batch, x,weights=varIdent(form=~1|Genotype),
                        na.action="na.omit", method="ML"))
        
        model_null <- do.call("lme", args=list(model_null.formula, x,
                        random=~1|Batch,weights=varIdent(form=~1|Genotype),
                        na.action="na.omit",  method="ML"))
        
        p.value <- (anova(model_genotype, model_null)$p[2])
    }else if(!keep_batch && !keep_equalvar){
        model_genotype <- do.call("gls", args = list(model_genotype.formula,
                        x,weights=varIdent(form=~1|Genotype),method="ML", na.action="na.omit"))
        
        model_null <- do.call("gls", args=list(model_null.formula,
                        x,weights=varIdent(form=~1|Genotype),method="ML", na.action="na.omit"))
        
        p.value <- (anova(model_genotype, model_null)$p[2])
    }else if(!keep_batch && keep_equalvar){
        model_genotype <- do.call("gls", args = list(model_genotype.formula,
                        x, method="ML", na.action="na.omit"))
        
        model_null <- do.call("gls", args=list(model_null.formula, x,
                        method="ML", na.action="na.omit"))
        
        p.value <- (anova(model_genotype, model_null)$p[2])
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
    
    ## Assign MM quality of fit
    result$model.output.quality <- testFinalModel(result)
    
    ## Parse modeloutput and choose output depending on model
    result$model.output.summary <- parserOutputSummary(result)
    
    return(result)
}
##------------------------------------------------------------------------------
## Parser model output summary and return in readable vector format
parserOutputSummary<-function(phenTestResult)
{
    result <- phenTestResult
    modeloutput_summary <- summary(result$model.output)
    genotype_estimate <- NA
    genotype_estimate_SE <- NA
    genotype_p_value <- NA
    
    gender_estimate <- NA
    gender_estimate_SE <- NA
    gender_p_value <- NA
    
    intercept_estimate <- NA
    intercept_estimate_SE <- NA
    weight_estimate <- NA
    weight_estimate_SE <- NA
    weight_p_value <- NA
    
    gender_FvKO_estimate <- NA
    gender_FvKO_SE <- NA
    gender_FvKO_p_value <- NA
    gender_MvKO_estimate <- NA
    gender_MvKO_SE <- NA
    gender_MvKO_p_value <- NA
    
    lengthoftable <- {
        table_length <- NA
        
        if (result$equation=="withWeight"){
            if(result$numberGenders==2){
                if((result$model.effect.gender
                                && result$model.effect.interaction)|
                        (!result$model.effect.gender && result$model.effect.interaction)){
                    table_length <- 5
                }else if(result$model.effect.gender &&
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
            if(result$numberGenders==2){
                if((result$model.effect.gender
                                && result$model.effect.interaction)|
                        (!result$model.effect.gender && result$model.effect.interaction)){
                    table_length <- 4
                }else if(!result$model.effect.gender &&
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
                if(result$model.effect.batch){
                    ## for mixed model
                    intercept_estimate = modeloutput_summary[["tTable"]][[1]]
                    intercept_estimate_SE = modeloutput_summary[["tTable"]][[(1+lengthoftable)]]
                    if((result$model.effect.gender && result$model.effect.interaction)
                            |( !result$model.effect.gender&& result$model.effect.interaction)){
                        gender_estimate=modeloutput_summary[["tTable"]][[2]]
                        gender_estimate_SE=modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                        gender_p_value= modeloutput_summary[["tTable"]][[(2+4*lengthoftable)]]
                        gender_FvKO_estimate= modeloutput_summary[["tTable"]][[3]]
                        gender_FvKO_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                        gender_FvKO_p_value=modeloutput_summary[["tTable"]][[(3+4*lengthoftable)]]
                        gender_MvKO_estimate=modeloutput_summary[["tTable"]][[4]]
                        gender_MvKO_SE=modeloutput_summary[["tTable"]][[(4+lengthoftable)]]
                        gender_MvKO_p_value=modeloutput_summary[["tTable"]][[(4+4*lengthoftable)]]
                    } else if( !result$model.effect.gender && !result$model.effect.interaction){
                        genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                        genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                        genotype_p_value =  modeloutput_summary[["tTable"]][[(2+4*lengthoftable)]]
                        
                    }else{
                        genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                        genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                        genotype_p_value =  modeloutput_summary[["tTable"]][[(2+4*lengthoftable)]]
                        gender_estimate=modeloutput_summary[["tTable"]][[3]]
                        gender_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                        gender_p_value= modeloutput_summary[["tTable"]][[(3+4*lengthoftable)]]
                    }
                    
                }else{
                    ## adaption for being a linear model rather than a mixed
                    ## model
                    intercept_estimate = modeloutput_summary[["tTable"]][[1]]
                    intercept_estimate_SE = modeloutput_summary[["tTable"]][[(1+lengthoftable)]]
                    
                    if((result$model.effect.gender && result$model.effect.interaction)
                            |( !result$model.effect.gender&& result$model.effect.interaction)){
                        
                        gender_estimate=modeloutput_summary[["tTable"]][[3]]
                        gender_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                        gender_p_value= modeloutput_summary[["tTable"]][[(3+3*lengthoftable)]]
                        gender_FvKO_estimate= modeloutput_summary[["tTable"]][[3]]
                        gender_FvKO_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                        gender_FvKO_p_value=modeloutput_summary[["tTable"]][[(3+3*lengthoftable)]]
                        gender_MvKO_estimate=modeloutput_summary[["tTable"]][[4]]
                        gender_MvKO_SE=modeloutput_summary[["tTable"]][[(4+lengthoftable)]]
                        gender_MvKO_p_value=modeloutput_summary[["tTable"]][[(4+3*lengthoftable)]]
                        
                        
                    } else if( !result$model.effect.gender && !result$model.effect.interaction){
                        
                        genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                        genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                        genotype_p_value =  modeloutput_summary[["tTable"]][[(2+3*lengthoftable)]]
                        
                    }else{
                        
                        genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                        genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                        genotype_p_value =  modeloutput_summary[["tTable"]][[(2+3*lengthoftable)]]
                        gender_estimate=modeloutput_summary[["tTable"]][[3]]
                        gender_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                        gender_p_value= modeloutput_summary[["tTable"]][[(3+3*lengthoftable)]]
                    }
                    
                    
                }
            },
            withWeight = {
                
                if(!result$model.effect.weight){
                    
                    ## If weight is not significant then the output is the
                    ## same as fitting model Eq1 and so no output is needed.
                    result$model.effect.batch=NA
                    
                }else{
                    
                    if(result$model.effect.batch){
                        
                        ## for mixed model
                        intercept_estimate = modeloutput_summary[["tTable"]][[1]]
                        intercept_estimate_SE = modeloutput_summary[["tTable"]][[(1+lengthoftable)]]
                        
                        if((result$model.effect.weight && result$model.effect.gender &&
                                        result$model.effect.interaction) |
                                (result$model.effect.weight &&
                                        !result$model.effect.gender &&
                                        result$model.effect.interaction)){
                            gender_estimate=modeloutput_summary[["tTable"]][[2]]
                            gender_estimate_SE=modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                            gender_p_value= modeloutput_summary[["tTable"]][[(2+4*lengthoftable)]]
                            gender_FvKO_estimate= modeloutput_summary[["tTable"]][[4]]
                            gender_FvKO_SE=modeloutput_summary[["tTable"]][[(4+lengthoftable)]]
                            gender_FvKO_p_value=modeloutput_summary[["tTable"]][[(4+4*lengthoftable)]]
                            gender_MvKO_estimate=modeloutput_summary[["tTable"]][[5]]
                            gender_MvKO_SE=modeloutput_summary[["tTable"]][[(5+lengthoftable)]]
                            gender_MvKO_p_value=modeloutput_summary[["tTable"]][[(5+4*lengthoftable)]]
                            weight_estimate=modeloutput_summary[["tTable"]][[3]]
                            weight_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                            weight_p_value=modeloutput_summary[["tTable"]][[(3+4*lengthoftable)]]
                            
                        } else if (result$model.effect.weight &&
                                !result$model.effect.gender &&
                                !result$model.effect.interaction){
                            
                            genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput_summary[["tTable"]][[(2+4*lengthoftable)]]
                            weight_estimate=modeloutput_summary[["tTable"]][[3]]
                            weight_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                            weight_p_value=modeloutput_summary[["tTable"]][[(3+4*lengthoftable)]]
                            
                        }else if (result$model.effect.weight &&
                                result$model.effect.gender &&
                                !result$model.effect.interaction){
                            
                            genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput_summary[["tTable"]][[(2+4*lengthoftable)]]
                            gender_estimate=modeloutput_summary[["tTable"]][[3]]
                            gender_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                            gender_p_value= modeloutput_summary[["tTable"]][[(3+4*lengthoftable)]]
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
                                        result$model.effect.gender &&
                                        result$model.effect.interaction )|
                                (result$model.effect.weight &&
                                        !result$model.effect.gender &&
                                        result$model.effect.interaction)){
                            gender_estimate=modeloutput_summary[["tTable"]][[2]]
                            gender_estimate_SE=modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                            gender_p_value= modeloutput_summary[["tTable"]][[(2+3*lengthoftable)]]
                            weight_estimate=modeloutput_summary[["tTable"]][[3]]
                            weight_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                            weight_p_value=modeloutput_summary[["tTable"]][[(3+3*lengthoftable)]]
                            gender_FvKO_estimate= modeloutput_summary[["tTable"]][[4]]
                            gender_FvKO_SE=modeloutput_summary[["tTable"]][[(4+lengthoftable)]]
                            gender_FvKO_p_value=modeloutput_summary[["tTable"]][[(4+3*lengthoftable)]]
                            gender_MvKO_estimate=modeloutput_summary[["tTable"]][[5]]
                            gender_MvKO_SE=modeloutput_summary[["tTable"]][[(5+lengthoftable)]]
                            gender_MvKO_p_value=modeloutput_summary[["tTable"]][[(5+3*lengthoftable)]]
                            
                        } else if (result$model.effect.weight &&
                                result$model.effect.gender &&
                                !result$model.effect.interaction){
                            
                            genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput_summary[["tTable"]][[(2+3*lengthoftable)]]
                            gender_estimate=modeloutput_summary[["tTable"]][[3]]
                            gender_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                            gender_p_value= modeloutput_summary[["tTable"]][[(3+3*lengthoftable)]]
                            weight_estimate=modeloutput_summary[["tTable"]][[4]]
                            weight_estimate_SE=modeloutput_summary[["tTable"]][[(4+lengthoftable)]]
                            weight_p_value=modeloutput_summary[["tTable"]][[(4+3*lengthoftable)]]
                            
                            
                        }else if (result$model.effect.weight &&
                                !result$model.effect.gender &&
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
            gender_estimate, gender_estimate_SE,  gender_p_value, 
            weight_estimate, weight_estimate_SE, weight_p_value, 
            intercept_estimate, intercept_estimate_SE, 
            gender_FvKO_estimate, gender_FvKO_SE, gender_FvKO_p_value,  
            gender_MvKO_estimate, gender_MvKO_SE, gender_MvKO_p_value)
    
    names(output) <- c("genotype_estimate", "genotype_estimate_SE", 
            "genotype_p_value", 
            "gender_estimate", "gender_estimate_SE", "gender_p_value", 
            "weight_estimate", "weight_estimate_SE", "weight_p_value", 
            "intercept_estimate", "intercept_estimate_SE", 
            "gender_FvKO_estimate", "gender_FvKO_SE", "gender_FvKO_p_value", 
            "gender_MvKO_estimate", "gender_MvKO_SE", "gender_MvKO_p_value")
    
    return(output)
}
##------------------------------------------------------------------------------