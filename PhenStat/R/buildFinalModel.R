# buildFinalModel.R contains buildFinalModel function

buildFinalModel <- function(object, result=NULL, equation=NULL, depVariable=NULL, pThreshold=0.05, keepList=NULL)

# Build final model based on the significance of different effects (see testDataset.R) 
# By default works with PhenTestResult object created by testDataset function.
# If someone would like to assign other TRUE/FALSE values to effects of the model then function 
# work with list of TRUE/FALSE values (keepList parameter).

{
    require(nlme)
    
    # Check PhenList object
    if(is(object,"PhenList")) {
        x <- object$phendata  
        
    } else {
        x <- as.data.frame(object)
    }
    
    # Check PhenTestResult object
    if(is(result,"PhenTestResult")) {
        if (is.null(depVariable)) depVariable <- result$depVariable
        if (is.null(equation)) equation <- result$equation
        keep_weight <- result$weightEffect
        keep_gender <- result$genderEffect
        keep_interaction <- result$interactionEffect
        keep_batch <- result$batchEffect
        keep_equalvar <- result$varianceEffect
    
        # Stop function if there are no enough needed input parameters
        if (is.null(depVariable)) stop("Please define dependant variable")
        if (is.null(equation)) stop("Please define equation: 'withWeight' or 'withoutWeight'")
        if (is.null(keep_batch) || is.null(keep_equalvar) || is.null(keep_gender) || is.null(keep_interaction)) 
        stop ("Please run function 'testDataset' first")
        if (result$equation!=equation) stop(paste("Tests have been done with another equation: ",result$equation))
        if (result$depVariable!=depVariable) stop(paste("Tests have been done for another dependant variable: ",result$depVariable))
    }
    else{
            # Stop function if there are no enough needed input parameters
            if (is.null(keepList) || length(keepList)!=5) 
                stop("Please define the values for 'keepList' list, where for each effect/part of the model TRUE/FALSE value defines to keep it in the model or not: 
                    'keepList=c(keepBatch,keepVariance,keepWeight,keepGender,keepInteraction)'")
            keep_weight <- keepList[3]
            keep_gender <- keepList[4]
            keep_interaction <- keepList[5]
            keep_batch <- keepList[1]
            keep_equalvar <- keepList[2]
        
            # Stop function if there are no enough needed input parameters
            if (is.null(depVariable)) stop("Please define dependant variable")
            if (is.null(equation)) stop("Please define equation: 'withWeight' or 'withoutWeight'")
        
            # Create start model
            model=buildStartModel(object,equation,depVariable,pThreshold,c(keep_batch,keep_equalvar))
        
            numberofgenders=length(levels(x$Gender))
        
            interactionTest=anova(model, type="marginal")$"p-value"[5]   
            
            # Create new PhenTestResult object using input parameters
            result <- new("PhenTestResult",list(modelOutput=model,depVariable=depVariable,equation=equation, 
                        batchEffect=keep_batch,varianceEffect=keep_equalvar,interactionEffect=keep_interaction,
                        interactionTestResult=interactionTest,genderEffect=keep_gender,weightEffect=keep_weight,
                        numberGenders=numberofgenders))
    }
   

    numberofgenders=result$numberGenders
    
    
    # Build final null model
    # Goal:  to test fixed effects of the model and based on the output build the final null model formula for later 
    # testing - as a null model it automatically excludes genotype and interaction term.
    # The anova function tests the fixed effects associated by treatment with a null hypothesis that the regression 
    # coefficients are equal to zero  and an alternative hypothesis that the regression coefficient are not equal to zero.
    # If the p-values of these tests are less than 0.05 we reject the null and accept the alternative that the are 
    # significant components of the model and should be included.
    # If no terms are significant a model can be build with just an intercept element this is specified as 
    # "model.formula <- as.formula(paste(depVariable, "~", "1"))"
    
    #Null model: genotype is not significant
    model_null.formula <- switch(equation,
            withWeight = {
                # Eq.2
                if(numberofgenders==2){
                    if(!keep_gender){
                        as.formula(paste(depVariable, "~", "Weight"))
                        
                    }else{
                        as.formula(paste(depVariable, "~", paste("Gender", "Weight", sep= "+")))
                    }
                }else{ 
                    as.formula(paste(depVariable, "~", "Weight"))
                }                 
            },
            withoutWeight = {
                # Eq.1
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
    
    #Alternative model: genotype is significant
    model_genotype.formula <- switch(equation,
            withWeight = {
                # Eq.2
                if(numberofgenders==2){
                    if ((keep_gender && keep_weight && keep_interaction)| (!keep_gender && keep_weight && keep_interaction)){
                        as.formula(paste(depVariable, "~", paste("Gender", "Genotype:Gender", "Weight", sep= "+")))
                        
                    } else if(keep_gender && keep_weight && !keep_interaction){
                        as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Weight", sep= "+")))
                        
                    } else if(!keep_gender && keep_weight && !keep_interaction){
                        as.formula(paste(depVariable, "~",paste("Genotype", "Weight", sep= "+")))
                    }
                    
                }else{
                    as.formula(paste(depVariable, "~", paste("Genotype", "Weight", sep="+")))
                }                
            },
            withoutWeight = {
                # Eq.1
                if(numberofgenders==2){
                    if (!keep_gender  && !keep_interaction){
                        as.formula(paste(depVariable, "~", "Genotype"))
                    } else if((keep_gender && keep_interaction)|(!keep_gender && keep_interaction)){
                        as.formula(paste(depVariable, "~", paste("Gender",  "Genotype:Gender", sep= "+")))
                    } else if(keep_gender && !keep_interaction){
                        as.formula(paste(depVariable, "~", paste("Genotype", "Gender", sep= "+")))
                    } 
                }else{
                    as.formula(paste(depVariable, "~", paste("Genotype")))     
                }                
            }            
    )
    
    
    # Test: genotype groups association with dependant variable
    # Null Hypothesis: genotypes are not associated with dependant variable 
    # Alternative Hypothesis: genotypes are associated with dependant variable 
    if(keep_batch && keep_equalvar){
        model_genotype=do.call("lme", args = list(model_genotype.formula, random=~1|Assay.Date, x, na.action="na.omit", method="ML"))
        model_null=do.call("lme", args=list(model_null.formula, x,random=~1|Assay.Date, na.action="na.omit",  method="ML"))
        p.value=(anova(model_genotype, model_null)$p[2])
    }else if(keep_batch && !keep_equalvar){
        model_genotype=do.call("lme", args = list(model_genotype.formula, random=~1|Assay.Date, x,weights=varIdent(form=~1|Genotype), na.action="na.omit", method="ML"))
        model_null=do.call("lme", args=list(model_null.formula, x, random=~1|Assay.Date,weights=varIdent(form=~1|Genotype), na.action="na.omit",  method="ML"))
        p.value=(anova(model_genotype, model_null)$p[2])
    }else if(!keep_batch && !keep_equalvar){
        model_genotype=do.call("gls", args = list(model_genotype.formula,  x,weights=varIdent(form=~1|Genotype),method="ML", na.action="na.omit"))
        model_null=do.call("gls", args=list(model_null.formula, x,weights=varIdent(form=~1|Genotype),method="ML", na.action="na.omit"))
        p.value=(anova(model_genotype, model_null)$p[2])
    }else if(!keep_batch && keep_equalvar){
        model_genotype=do.call("gls", args = list(model_genotype.formula,  x, method="ML", na.action="na.omit"))
        model_null=do.call("gls", args=list(model_null.formula, x,method="ML", na.action="na.omit"))
        p.value=(anova(model_genotype, model_null)$p[2])
    }
    
    # Final model version with na.exclude and REML method
    if(keep_batch && keep_equalvar){
        model_genotype=do.call("lme", args = list(model_genotype.formula, random=~1|Assay.Date, x, na.action="na.exclude", method="REML"))
    }else if(keep_batch && !keep_equalvar){
        model_genotype=do.call("lme", args = list(model_genotype.formula, random=~1|Assay.Date, x,weights=varIdent(form=~1|Genotype), na.action="na.exclude", method="REML"))
    }else if(!keep_batch && !keep_equalvar){
        model_genotype=do.call("gls", args = list(model_genotype.formula,  x,weights=varIdent(form=~1|Genotype), na.action="na.exclude"))
    }else if(!keep_batch && keep_equalvar){
        model_genotype=do.call("gls", args = list(model_genotype.formula,  x, na.action="na.exclude"))
    }

   
    # Store the results      
    result$modelOutput=model_genotype
    result$genotypeEffect=p.value
    result$modelFormula.null=model_null.formula
    result$modelFormula.genotype=model_genotype.formula
    
    # Assign MM quality of fit
    MM_fitquality=Diagnostictest(x,result)
    
    result$MM_fitquality=MM_fitquality

    return(result)
    
    
    
}