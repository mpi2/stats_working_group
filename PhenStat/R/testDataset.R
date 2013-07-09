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

# testDataset.R contains testDataset, buildStartModel and modelFormula functions

testDataset <- function(object, equation="withWeight", depVariable=NULL, pThreshold=0.05, method="MM")

# Create start model and modify it after testing of different hypothesis.
# TRUE/FALSE values assigned to effects of the model are stored in the PhenTestResult object for the further final model build.

# The testable effects are: batch effect (random effects significance), variance effect (TRUE if residual variances for 
# genotype groups are homogeneous and FALSE if they are heterogeneous), interaction effect (genotype by gender 
# interaction significance) plus interaction test anova results, gender effect (gender significance), weigth effect 
# (weigth significance).

# TODO
# check dataset for dep. variable


{
    require(nlme)
    
    # Check PhenList object
    if(is(object,"PhenList")) {
        x <- object$dataset    
        
    } else {
        # generic message about dataset format - advise...
        x <- as.data.frame(object)
        x <- checkDataset(x)
    }
    
    
    # Test: depVariable is continuous variable
    columnOfInterest <- x[,c(depVariable)]
    if(is.numeric(columnOfInterest)){
        if (length(unique(columnOfInterest))/length(columnOfInterest)<0.05) 
            message(paste("Warning: Dependent variable '",depVariable,"' is numeric but seemed to be cathegorical. Fisher exact test can be better way to do the analysis than mixed models.",sep="")) 
    }
    else stop(paste("Error: Dependent variable '",depVariable,"' is not numeric. Please run Fisher exact test for the analysis of this dependent variable.",sep="")) 
    
    
    # Mixed Models
    if (method=="MM"){ 


        # Stop function if there are no enough needed input parameters
        if (is.null(depVariable)) stop("Please define dependant variable")
        
        numberofgenders=length(levels(x$Gender))
        
        # Start model formula: homogenous residual variance, genotype and sex interaction included  
        model.formula = modelFormula(equation,numberofgenders, depVariable)
        
        if ('Batch' %in% colnames(dataset)){
        
            # MM fit of model formula (with random effects)
            # Model 1 
            model_MM = do.call("lme", args = list(model.formula, random=~1|Batch, data = x, na.action="na.omit", method="REML"))
            # GLS fit of model formula  (no random effects)
            # Model 1A
            model_withoutbatch <- do.call("gls", args=list(model.formula, x, na.action="na.omit"))
            # Test: the random effects associated with batch intercepts can be ommited from model
            # Hypothesis 1
            # Null Hypothesis: variance of batch = 0
            # Alternative Hypothesis: variance of batch > 0 
            p.value.batch <-(anova(model_MM, model_withoutbatch)$p[2])/2
            # The result of the test for Hypothesis 1 will help to select the structure for random effects
            keep_batch= p.value.batch<pThreshold
            
            # MM fit of model formula with heterogeneous residual variances for genotype groups
            # Model 1 assumes homogeneous residual variances
            # Model 2 with heterogeneous residual variances 
            model_hetvariance= do.call("lme", args=list(model.formula, random=~1|Batch, x, weights=varIdent(form=~1|Genotype), na.action="na.omit", method="REML"))
            # Test: the variance of the residuals is the same (homogeneous) for all genotype groups
            # Hypothesis 2
            # Null Hypothesis: all residual variances are equal
            # Alternative Hypothesis: at least one pair of residual variances is not equal    
            p.value.variance=(anova(model_MM, model_hetvariance)$p[2])
            # The result of the test for Hypothesis 2 will help to select a covariance structure for the residuals
            keep_equalvar= p.value.variance>pThreshold
        }
        else {
            # No Batch effects
            keep_batch = FALSE
            keep_equalvar = FALSE
        }
        
        
        # Model fit is selected according to test results    
        if(keep_batch && keep_equalvar){
            # Model 1
            model= model_MM
        }else if(keep_batch && !keep_equalvar){
            # Model 2
            model= model_hetvariance
        }else if(!keep_batch && keep_equalvar){
            # Model 1A
            model= model_withoutbatch
        }else if(!keep_batch && !keep_equalvar){
            # Modify model 1A to heterogeneous residual variances
            model= do.call("gls", args=list(model.formula, weights=varIdent(form=~1|Genotype), x, na.action="na.omit"))
        }
    

        # Tests for significance of fixed effects using TypeI F-test from anova functionality by using selected model
        anova_results = anova(model, type="marginal")$"p-value" < pThreshold
        if(numberofgenders==2){
            # Result of the test for gender significance (fixed effect 1.)
            keep_gender = anova_results[3]
            # Eq.2
            if (equation=="withWeight"){   
                # Result of the test for weight significance  (fixed effect 3.)        
                keep_weight = anova_results[4]
                # Result of the test for genotype by gender interaction significance (fixed effect 2.)
                keep_interaction = anova_results[5]
                
                # Technical results needed for the output
                # Interaction test results are kept for the output 
                interactionTest=anova(model, type="marginal")$"p-value"[5]           
            }
            # Eq.1
            else{
                # Result of the test for weight significance  (fixed effect 3.) 
                # It's FALSE since here the equasion 1 is used - without weight effect
                keep_weight = FALSE
                # Result of the test for genotype by gender interaction significance (fixed effect 2.)
                keep_interaction = anova_results[4]
                # Interaction test results are kept for the output
                interactionTest=anova(model, type="marginal")$"p-value"[4]      
            } 
        }
        
        result <- new("PhenTestResult",list(model.output=model,depVariable=depVariable,equation=equation, 
                        model.effect.batch=keep_batch,model.effect.variance=keep_equalvar,model.effect.interaction=keep_interaction,
                        model.output.interaction=interactionTest,model.effect.gender=keep_gender,model.effect.weight=keep_weight,
                        numberGenders=numberofgenders))
        #pThreshold should be stored
    }
    else
    # Fisher Exact Test placeholder
    {
        message("Fisher Exact Test")
        result <- NULL
    }
    
    return(result)   
}

buildStartModel <- function(object, equation="withWeight", depVariable, pThreshold=0.05, keepList)
# If someone would like to assign other TRUE/FALSE values to effects of the model then start model is build by using
# buildStartModel function  
{
    require(nlme)
    
    # Check object
    if(is(object,"PhenList")) {
        x <- object$dataset    
        
    } else {
        x <- as.data.frame(object)
    }
    
    numberofgenders=length(levels(x$Gender))
    
    # User's values for effects
    keep_batch <- keepList[1]
    keep_equalvar <- keepList[2]
    
    # Start model formula: homogenous residual variance, genotype and sex interaction included      
    model.formula = modelFormula(equation, numberofgenders, depVariable)
    
    # MM fit of model formula (with random effects)
    # Model 1 
    model_MM = do.call("lme", args = list(model.formula, random=~1|Batch, data = x, na.action="na.omit", method="REML"))
    # GLS fit of model formula  (no random effects)
    # Model 1A
    model_withoutbatch <- do.call("gls", args=list(model.formula, x, na.action="na.omit"))
    # MM fit of model formula with heterogeneous residual variances for genotype groups
    # Model 1 assumes homogeneous residual variances
    # Model 2 with heterogeneous residual variances 
    model_hetvariance= do.call("lme", args=list(model.formula, random=~1|Batch, x, weights=varIdent(form=~1|Genotype), 
                    na.action="na.omit", method="REML"))
    
    # Model fit is selected according to test results    
    if(keep_batch && keep_equalvar){
        # Model 1
        model= model_MM
    }else if(keep_batch && !keep_equalvar){
        # Model 2
        model= model_hetvariance
    }else if(!keep_batch && keep_equalvar){
        # Model 1A
        model= model_withoutbatch
    }else if(!keep_batch && !keep_equalvar){
        # Modify model 1A to heterogeneous residual variances
        model= do.call("gls", args=list(model.formula, weights=varIdent(form=~1|Genotype), x, na.action="na.omit"))
    }
    
    return(model)
}


modelFormula <- function(equation, numberofgenders, depVariable)
# Creats formula for the start model based on equation and number of genders in the data
{
    
    model.formula <- switch(equation,
            # Eq.2
            withWeight = {
                # Fixed effects: 1) Gender 2) Genotype by Gender interaction 3) Weight
                if(numberofgenders==2){                            
                    as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Genotype*Gender","Weight", sep= "+")))
                }else{ 
                    as.formula(paste(depVariable, "~", paste("Genotype", "Weight",  sep= "+"))) 
                }         
            },
            # Eq.1 
            withoutWeight = {
                # Fixed effects: 1) Gender 2) Genotype by Gender interaction
                if(numberofgenders==2){
                    as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Genotype*Gender", sep= "+")))
                }else{ 
                    as.formula(paste(depVariable, "~", paste("Genotype",  sep= "+"))) 
                } 
            }
            )
    return(model.formula)
}