testDataset <- function(object, equation="withWeight", depVariable=NULL, pThreshold=0.05)
# Function to test batch, variance, gender and interaction effects
# A full loaded model formula (ie one with all variables of interest) is assembled.  Weight cannot be in the model as a dependent variable if it is the independent variable - hence the if else questions to build the model.
# Then the model formula are used to build two models one with batch included as a random effect and hence uses a mixed model and one where batch is not included and uses linear regression
# The two models are compared via anova to test the null hypothesis that batch is not significant (Ho) where the alternative hypothesis (Ha) is that batch is significant.  If the p value <0.05,  then the we reject the null hypothesis and accept the alternative and say we need to keep batch
# the p-value for batch test is divided by 2 as the test is on the boundary of the parameter space for a variance and the null distribution of the likelihood ratio test statistics follows a mixture of chi-squared distributions with equal weight of 0.5.
# As the output is yes or no on keep_batch then if p<0.05,  we say batch is significant  and a mixed model approach should be used in subsequent work.

{
    require(nlme)
    #    Check object
    if(is(object,"PhenList")) {
        x <- object$phendata
        #if (is.null(depVariable)) depVariable <- object$depVariable
        #if (is.null(equation)) equation <- object$equation
        
    } else {
        x <- as.data.frame(object)
    }
            
    if (is.null(depVariable)) stop("Please define dependant variable")
    if (is.null(equation)) stop("Please define equation: 'withWeight' or 'withoutWeight'")
    
    #    Check method    
    #method <- match.arg(method)
    numberofgenders=length(levels(x$Gender))
            
    model.formula <- switch(equation,
                    withWeight = {
                        
                        if(numberofgenders==2){
                            as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Genotype*Gender","Weight", sep= "+")))
                        }else{ 
                            as.formula(paste(depVariable, "~", paste("Genotype", "Weight",  sep= "+"))) 
                        } 
                        
                    },
                    withoutWeight = {
                        
                        if(numberofgenders==2){
                            as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Genotype*Gender", sep= "+")))
                        }else{ 
                            as.formula(paste(depVariable, "~", paste("Genotype",  sep= "+"))) 
                        } 
                        
                    }
                    
    )

    #Testing batch
    model_MM = do.call("lme", args = list(model.formula, random=~1|Assay.Date, data = x, na.action="na.omit", method="REML"))
    model_withoutbatch <- do.call("gls", args=list(model.formula, x, na.action="na.omit"))
    p.value.batch <-(anova(model_MM, model_withoutbatch)$p[2])/2
    keep_batch= p.value.batch<pThreshold
    
    #Testing variance
    model_hetvariance= do.call("lme", args=list(model.formula, random=~1|Assay.Date, x, weights=varIdent(form=~1|Genotype), na.action="na.omit", method="REML"))
    p.value.variance=(anova(model_MM, model_hetvariance)$p[2])
    keep_equalvar= p.value.variance>pThreshold
    
    # Build model following variance and batch test for testing of fixed effects
    # The keep_batch and keep_equal variance functions are called to build the correct models 
    # via a series of if else rules.
    
    if(keep_batch && keep_equalvar){
        model= model_MM
    }else if(keep_batch && !keep_equalvar){
        model= model_hetvariance
    }else if(!keep_batch && keep_equalvar){
        model= model_withoutbatch
    }else if(!keep_batch && !keep_equalvar){
        model= do.call("gls", args=list(model.formula, weights=varIdent(form=~1|Genotype), x, na.action="na.omit"))
    }

    
    #Testing gender, interaction
    anova_results = anova(model, type="marginal")$"p-value" < 0.05
    

    if(numberofgenders==2){
        keep_gender = anova_results[3]
        if (equation=="withWeight"){            
            keep_weight = anova_results[4]
            keep_interaction = anova_results[5]
            interactionTest=anova(model, type="marginal")$"p-value"[5]
            
            # table length
            if(numberofgenders==2){
                if((keep_gender && keep_interaction)|(!keep_gender && keep_interaction)){
                    table_length=5
                }else if(keep_gender && !keep_interaction){
                    table_length=4
                }else {
                    table_length=3
                }  
            }else{
                table_length=3
            }

        }
        else{
            keep_weight = FALSE
            keep_interaction = anova_results[4]
            interactionTest=anova(model, type="marginal")$"p-value"[4]
            
            #table length 
            #How many elements will be on the table based on how the model used.  
            #However the model used is not a direct additive system of the yes calls (ie you can't have an interaction if you don't have gender in the model).
            if(numberofgenders==2){
                if((keep_gender && keep_interaction)|(!keep_gender && keep_interaction)){
                    table_length=4
                }else if(!keep_gender && !keep_interaction){
                    table_length=2
                }else{
                    table_length=3
                }  
            }else{
                table_length=2
            }


        } 
    }

    
    
    #    Output


        result <- new("PhenTestResult",list(modelOutput=NULL,depVariable=depVariable,equation=equation, 
                        batchEffect=keep_batch,varianceEffect=keep_equalvar,interactionEffect=keep_interaction,
                        interactionTest=interactionTest,genderEffect=keep_gender,weightEffect=keep_weight,outputLength=table_length))

        
        return(result)
   
}