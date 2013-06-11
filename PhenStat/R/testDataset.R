testDataset <- function(object, equation="withWeight", depVariable=NULL, pThreshold=0.05)

# Start model is created and modified by using different hypothesis: random effects significance etc. 
# Significance is assigned to fixed effects  
# As the result important effects are stored in the PhenTestResult object for the further final model build.
# These effects are: batch effect (random effects significance), variance effect (TRUE if residual variances for 
# genotype groups are homogeneous and FALSE if they are heterogeneous), interaction effect (genotype by gender 
# interaction significance) plus interaction test anova results, gender effect (gender significance), weigth effect 
# (weigth significance) and finally, output length estimation based on tests results)
# The idea here is regardless the test results, all effects can be change if needed before the final model build

{
    require(nlme)

    #    Check object
    if(is(object,"PhenList")) {
        x <- object$phendata    
        
    } else {
        x <- as.data.frame(object)
    }
           
    # Stop function if there are no enough needed input parameters
    if (is.null(depVariable)) stop("Please define dependant variable")
    if (is.null(equation)) stop("Please define equation: 'withWeight' or 'withoutWeight'")
    

    numberofgenders=length(levels(x$Gender))
       
    # Start model formula: homogenous residual variance, genotype and sex interaction included      
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

    # MM fit of model formula (with random effects)
    # Model 1 
    model_MM = do.call("lme", args = list(model.formula, random=~1|Assay.Date, data = x, na.action="na.omit", method="REML"))
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
    model_hetvariance= do.call("lme", args=list(model.formula, random=~1|Assay.Date, x, weights=varIdent(form=~1|Genotype), na.action="na.omit", method="REML"))
    # Test: the variance of the residuals is the same (homogeneous) for all genotype groups
    # Hypothesis 2
    # Null Hypothesis: all residual variances are equal
    # Alternative Hypothesis: at least one pair of residual variances is not equal    
    p.value.variance=(anova(model_MM, model_hetvariance)$p[2])
    # The result of the test for Hypothesis 2 will help to select a covariance structure for the residuals
    keep_equalvar= p.value.variance>pThreshold
    

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

    anova_results = anova(model, type="marginal")$"p-value" < 0.05
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

        result <- new("PhenTestResult",list(modelOutput=model,depVariable=depVariable,equation=equation, 
                        batchEffect=keep_batch,varianceEffect=keep_equalvar,interactionEffect=keep_interaction,
                        interactionTestResult=interactionTest,genderEffect=keep_gender,weightEffect=keep_weight,
                        numberGenders=numberofgenders))
        
        return(result)
   
}

buildStartModel <- function(object, equation="withWeight", depVariable, pThreshold=0.05, keepList)
{
    require(nlme)
    
    #    Check object
    if(is(object,"PhenList")) {
    x <- object$phendata    
    
    } else {
    x <- as.data.frame(object)
    }
    
    numberofgenders=length(levels(x$Gender))

    keep_batch <- keepList[1]
    keep_equalvar <- keepList[2]

    # Start model formula: homogenous residual variance, genotype and sex interaction included      
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

    # MM fit of model formula (with random effects)
    # Model 1 
    model_MM = do.call("lme", args = list(model.formula, random=~1|Assay.Date, data = x, na.action="na.omit", method="REML"))
    # GLS fit of model formula  (no random effects)
    # Model 1A
    model_withoutbatch <- do.call("gls", args=list(model.formula, x, na.action="na.omit"))
    # MM fit of model formula with heterogeneous residual variances for genotype groups
    # Model 1 assumes homogeneous residual variances
    # Model 2 with heterogeneous residual variances 
    model_hetvariance= do.call("lme", args=list(model.formula, random=~1|Assay.Date, x, weights=varIdent(form=~1|Genotype), 
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
