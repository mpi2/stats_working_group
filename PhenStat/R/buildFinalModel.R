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
        x <- object$dataset  
        
    } else {
        x <- as.data.frame(object)
        x <- checkDataset(x)
    }
    
    # Check PhenTestResult object
    if(is(result,"PhenTestResult")) {
        if (is.null(depVariable)) depVariable <- result$depVariable
        if (is.null(equation)) equation <- result$equation
        keep_weight <- result$model.effect.weight
        keep_gender <- result$model.effect.gender
        keep_interaction <- result$model.effect.interaction
        keep_batch <- result$model.effect.batch
        keep_equalvar <- result$model.effect.variance
    
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
            result <- new("PhenTestResult",list(model.output=model,depVariable=depVariable,equation=equation, 
                        model.effect.batch=keep_batch,model.effect.variance=keep_equalvar,model.effect.interaction=keep_interaction,
                        model.output.interaction=interactionTest,model.effect.gender=keep_gender,model.effect.weight=keep_weight,
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
        model_genotype=do.call("lme", args = list(model_genotype.formula, random=~1|Batch, x, na.action="na.omit", method="ML"))
        model_null=do.call("lme", args=list(model_null.formula, x,random=~1|Batch, na.action="na.omit",  method="ML"))
        p.value=(anova(model_genotype, model_null)$p[2])
    }else if(keep_batch && !keep_equalvar){
        model_genotype=do.call("lme", args = list(model_genotype.formula, random=~1|Batch, x,weights=varIdent(form=~1|Genotype), na.action="na.omit", method="ML"))
        model_null=do.call("lme", args=list(model_null.formula, x, random=~1|Batch,weights=varIdent(form=~1|Genotype), na.action="na.omit",  method="ML"))
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
        model_genotype=do.call("lme", args = list(model_genotype.formula, random=~1|Batch, x, na.action="na.exclude", method="REML"))
    }else if(keep_batch && !keep_equalvar){
        model_genotype=do.call("lme", args = list(model_genotype.formula, random=~1|Batch, x,weights=varIdent(form=~1|Genotype), na.action="na.exclude", method="REML"))
    }else if(!keep_batch && !keep_equalvar){
        model_genotype=do.call("gls", args = list(model_genotype.formula,  x,weights=varIdent(form=~1|Genotype), na.action="na.exclude"))
    }else if(!keep_batch && keep_equalvar){
        model_genotype=do.call("gls", args = list(model_genotype.formula,  x, na.action="na.exclude"))
    }

   
    # Store the results      
    result$model.output=model_genotype
    result$model.output.genotype.nulltest.pVal=p.value
    result$model.formula.null=model_null.formula
    result$model.formula.genotype=model_genotype.formula
    
    # Assign MM quality of fit
    result$model.output.quality=testFinalModel(object,result)
    
    # Parse modeloutput and choose output depending on model 
    result$model.output.summary = parserOutputSummary(result)

    return(result)  
    
}

parserOutputSummary<-function(result)

{
    modeloutput_summary = summary(result$model.output)
    genotype_estimate =NA
    genotype_estimate_SE =NA
    genotype_p_value =NA
    
    gender_estimate=NA
    gender_estimate_SE=NA
    gender_p_value=NA
    
    intercept_estimate = NA
    intercept_estimate_SE =NA
    weight_estimate=NA
    weight_estimate_SE=NA
    weight_p_value=NA
    
    gender_FvKO_estimate=NA
    gender_FvKO_SE=NA
    gender_FvKO_p_value=NA
    gender_MvKO_estimate=NA
    gender_MvKO_SE=NA
    gender_MvKO_p_value=NA
    
    lengthoftable=outputLength(result)
    
    switch(result$equation,
            withoutWeight = {
                if(result$model.effect.batch){                   
                    #for mixed model 
                    intercept_estimate = modeloutput_summary[["tTable"]][[1]]
                    intercept_estimate_SE = modeloutput_summary[["tTable"]][[(1+lengthoftable)]]
                    if((result$model.effect.gender && result$model.effect.interaction) |( !result$model.effect.gender&& result$model.effect.interaction)){
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
                    #adaption for being a linear model rather than a mixed model
                    intercept_estimate = modeloutput_summary[["tTable"]][[1]]
                    intercept_estimate_SE = modeloutput_summary[["tTable"]][[(1+lengthoftable)]]
                    
                    if((result$model.effect.gender && result$model.effect.interaction) |( !result$model.effect.gender&& result$model.effect.interaction)){
                        
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
                    
                    #If weight is not significant then the output is the same as fitting model Eq1 and so no output is needed. 
                    result$model.effect.batch=NA
                    
                }else{
                    
                    if(result$model.effect.batch){
                        
                        #for mixed model 
                        intercept_estimate = modeloutput_summary[["tTable"]][[1]]
                        intercept_estimate_SE = modeloutput_summary[["tTable"]][[(1+lengthoftable)]]
                        
                        if((result$model.effect.weight && result$model.effect.gender && result$model.effect.interaction) | (result$model.effect.weight && !result$model.effect.gender&& result$model.effect.interaction)){
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
                            
                        } else if (result$model.effect.weight && !result$model.effect.gender && !result$model.effect.interaction){    
                            
                            genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput_summary[["tTable"]][[(2+4*lengthoftable)]]
                            weight_estimate=modeloutput_summary[["tTable"]][[3]]
                            weight_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                            weight_p_value=modeloutput_summary[["tTable"]][[(3+4*lengthoftable)]]
                            
                        }else if (result$model.effect.weight && result$model.effect.gender && !result$model.effect.interaction){
                            
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
                        #adaption for being a linear model rather than a mixed model
                        intercept_estimate = modeloutput_summary[["tTable"]][[1]]
                        intercept_estimate_SE = modeloutput_summary[["tTable"]][[(1+lengthoftable)]]
                        
                        if((result$model.effect.weight && result$model.effect.gender && result$model.effect.interaction )|(result$model.effect.weight && !result$model.effect.gender&& result$model.effect.interaction)){
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
                            
                        } else if (result$model.effect.weight && result$model.effect.gender && !result$model.effect.interaction){
                            
                            genotype_estimate = modeloutput_summary[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput_summary[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput_summary[["tTable"]][[(2+3*lengthoftable)]]
                            gender_estimate=modeloutput_summary[["tTable"]][[3]]
                            gender_estimate_SE=modeloutput_summary[["tTable"]][[(3+lengthoftable)]]
                            gender_p_value= modeloutput_summary[["tTable"]][[(3+3*lengthoftable)]]
                            weight_estimate=modeloutput_summary[["tTable"]][[4]]
                            weight_estimate_SE=modeloutput_summary[["tTable"]][[(4+lengthoftable)]]
                            weight_p_value=modeloutput_summary[["tTable"]][[(4+3*lengthoftable)]]
                            
                            
                        }else if (result$model.effect.weight && !result$model.effect.gender && !result$model.effect.interaction){    
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
    output = c(genotype_estimate, genotype_estimate_SE,  genotype_p_value,gender_estimate, gender_estimate_SE,  gender_p_value, weight_estimate, weight_estimate_SE, weight_p_value, intercept_estimate, intercept_estimate_SE, gender_FvKO_estimate, gender_FvKO_SE,   gender_FvKO_p_value,  gender_MvKO_estimate,  gender_MvKO_SE,  gender_MvKO_p_value)
    names = c("genotype_estimate", "genotype_estimate_SE", "genotype_p_value", "gender_estimate", "gender_estimate_SE", "gender_p_value", "weight_estimate", 
            "weight_estimate_SE", "weight_p_value", "intercept_estimate", "intercept_estimate_SE", "gender_FvKO_estimate", "gender_FvKO_SE", 
            "gender_FvKO_p_value", "gender_MvKO_estimate", "gender_MvKO_SE", "gender_MvKO_p_value")
    names(output) = names
    return(output)
}