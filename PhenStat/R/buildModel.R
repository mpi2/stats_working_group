buildModel <- function(object, result=NULL, equation=NULL, depVariable=NULL, pThreshold=0.05, keep_list=NULL)
#Function testing_batch: function to test batch
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
        if (!is(x$Genotype)) stop("Genotype values are not defined")   
        
    } else {
        x <- as.data.frame(object)
        if (!is(x$Genotype)) stop("Genotype values are not defined")   
    }
    
    if(is(result,"PhenTestResult")) {
        if (is.null(depVariable)) depVariable <- result$depVariable
        if (is.null(equation)) equation <- result$equation
        keep_weight <- result$weightEffect
        keep_gender <- result$genderEffect
        keep_interaction <- result$interactionEffect
        keep_batch <- result$batchEffect
        keep_equalvar <- result$varianceEffect
    
        if (is.null(depVariable)) stop("Please define dependant variable")
        if (is.null(equation)) stop("Please define equation: 'withWeight' or 'withoutWeight'")
        if (is.null(keep_batch) || is.null(keep_equalvar) || is.null(keep_gender) || is.null(keep_interaction)) 
        stop ("Please run function 'testData' first")
        if (result$equation!=equation) stop(paste("Tests have been done with another equation: ",result$equation))
        if (result$depVariable!=depVariable) stop(paste("Tests have been done for another dependant variable: ",result$depVariable))
    }
    else{
            if (is.null(keep_list) || length(keep_list)!=5) stop("Please define the values for tests: 
                    'keep_list=c(keep_batch,keep_equalvar,keep_weight,keep_gender,keep_interaction)'")
            keep_weight <- keep_list[3]
            keep_gender <- keep_list[4]
            keep_interaction <- keep_list[5]
            keep_batch <- keep_list[1]
            keep_equalvar <- keep_list[2]
            result <- new("PhenTestResult",list(modelOutput=NULL,depVariable=depVariable,equation=equation, 
                        batchEffect=keep_batch,varianceEffect=keep_equalvar,interactionEffect=keep_interaction,
                        genderEffect=keep_gender,weightEffect=keep_weight))
    }
   

    numberofgenders=length(levels(x$Gender))
    
    
    # Build final null model
    # Goal:  to test fixed effects of the model and based on the output build the final null model formula for later 
    # testing - as a null model it automatically excludes genotype and interaction term.
    # The anova function tests the fixed effects associated by treatment with a null hypothesis that the regression 
    # coefficients are equal to zero  and an alternative hypothesis that the regression coefficient are not equal to zero.
    # If the p-values of these tests are less than 0.05 we reject the null and accept the alternative that the are 
    # significant components of the model and should be included.
    # If no terms are significant a model can be build with just an intercept element this is specified as 
    # "model.formula <- as.formula(paste(depVariable, "~", "1"))"
    
    #Null model
    model_null.formula <- switch(equation,
            withWeight = {
                
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
    
    #Genotype model
    model_genotype.formula <- switch(equation,
            withWeight = {
                
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
    
    #finalmodel version with na.exclude

    if(keep_batch && keep_equalvar){
        model_genotype=do.call("lme", args = list(model_genotype.formula, random=~1|Assay.Date, x, na.action="na.exclude", method="REML"))
    }else if(keep_batch && !keep_equalvar){
        model_genotype=do.call("lme", args = list(model_genotype.formula, random=~1|Assay.Date, x,weights=varIdent(form=~1|Genotype), na.action="na.exclude", method="REML"))
    }else if(!keep_batch && !keep_equalvar){
        model_genotype=do.call("gls", args = list(model_genotype.formula,  x,weights=varIdent(form=~1|Genotype), na.action="na.exclude"))
    }else if(!keep_batch && keep_equalvar){
        model_genotype=do.call("gls", args = list(model_genotype.formula,  x, na.action="na.exclude"))
    }

    # MM fit quality
    a=levels(x$Genotype)
    numberofgenders=length(levels(x$Gender))
    #withWeight and weight is not significant
    if(!keep_weight && equation=="withWeight"){
        testresults=c(a[1], NA, a[2], NA, NA, NA)
        
    }else{    
        #withoutWeight or weight significant
        require(nortest)
        res=resid(model_genotype)
        data_all= data.frame(x, res)
        genotype_no=length(a)
        data_all[, c("Gender", "Assay.Date")] = lapply(data_all[, c("Gender", "Assay.Date")], factor)
        No_batches=nlevels(data_all$Assay.Date)
        outputnumeric=is.numeric(model_genotype$apVar)
        
        if(keep_batch && No_batches >7 && outputnumeric){
            blups=ranef(model_genotype)
            blups_test= cvm.test(blups [ ,1])$p.value
            sdests = exp(attr(model_genotype$apVar, "Pars"))           #extract variance estimates
            Zbat = model.matrix(~ Assay.Date, model.frame( ~ Assay.Date, model_genotype$groups))    #create random effects design matrix
            ycov = (Zbat %*% t(Zbat)) * sdests["reStruct.Assay.Date"]^2 + diag(rep(1,nrow(model_genotype$groups))) * sdests["lSigma"]^2    #create estimated cov(y)
            Lt = chol(solve(ycov))  #Cholesky decomposition of inverse of cov(y) (see Houseman '04 eq. (2))
            rotres = Lt %*%  model_genotype$residuals[, "fixed"]    #rotated residuals
            rotated_residual_test=cvm.test(rotres)$p.value
        }else{
            blups_test=NA
            rotated_residual_test=NA
        }   
        
        Gp1 = subset(data_all, data_all$Genotype==a[1])
        Gp2 = subset(data_all, data_all$Genotype==a[2])
        No_Gp1 = sum(is.finite(Gp1[ , depVariable])) 
        No_Gp2 = sum(is.finite(Gp2[ , depVariable])) 
        
        if(No_Gp1>7){
            gp1_norm_res= cvm.test(Gp1$res)$p.value
        }else{
            gp1_norm_res= NA
        }    
        
        if(No_Gp2>7){
            gp2_norm_res= cvm.test(Gp2$res)$p.value
        }else{
            gp2_norm_res= NA
        }    
        
        testresults=c(a[1], gp1_norm_res, a[2], gp2_norm_res, blups_test, rotated_residual_test)
        
    }

    message("Null model formula:")
    message(model_null.formula)
    message("Genotype model formula:")
    message(model_genotype.formula)
    message("Genotype effect:")
    message(p.value)
        
    result$modelOutput=model_genotype
    result$genotypeEffect=p.value
    result$model.null=model_null.formula
    result$model.genotype=model_genotype.formula
    result$MM_fitquality=testresults
    #Final Model
    #message("Model formula from final model function:")
    return(result)
    
    
    
}