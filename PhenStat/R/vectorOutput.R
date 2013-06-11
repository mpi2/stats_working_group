# vectorOutput.R contains vectorOutput and outputLength functions

vectorOutput <- function(result)
# Wrapper to prepare the output of the modeling and testing results in vector form
{
    if(!is(result,"PhenTestResult")) {
        stop ("Please provide PhenTestResult object as an argument")
    }
    
    if (is.null(result$modelFormula.genotype)) stop("There are no results to wrap. Please run function 'buildFinalModel' first")

    MM_fitquality=result$MM_fitquality
    depVariable=result$depVariable
    
    # Fitted final genotype model
    modeloutput=result$modelOutput 
    
    keep_batch=result$batchEffect
    variance_test=result$varianceEffect
    Nulltest_genotype_pvalue=result$genotypeEffect

    Gender_sig=result$genderEffect
    Weight_sig=result$weightEffect
    Interaction_sig=result$interactionEffect
    
    Interaction_test=result$interactionTestResult
    
    # Problem depending on the length of the table where we grab values. 
    # There is one less column when batch is not significant.
    lengthoftable=outputLength(result)
    
 
    # The table is organised as a vector of values,  
    # Ordered going down each column in the summary table.
    values <- switch(result$equation,
            withoutWeight = {
                if(keep_batch){
                    #for mixed model 
                    intercept_estimate = modeloutput[["tTable"]][[1]]
                    intercept_estimate_SE = modeloutput[["tTable"]][[(1+lengthoftable)]]
                    weight_estimate=NA
                    weight_estimate_SE=NA
                    weight_p_value=NA
                    if((Gender_sig && Interaction_sig) |( !Gender_sig&& Interaction_sig)){
                        genotype_estimate =NA
                        genotype_estimate_SE =NA
                        genotype_p_value =NA
                        gender_estimate=modeloutput[["tTable"]][[2]]
                        gender_estimate_SE=modeloutput[["tTable"]][[(2+lengthoftable)]]
                        gender_p_value= modeloutput[["tTable"]][[(2+4*lengthoftable)]]
                        gender_FvKO_estimate= modeloutput[["tTable"]][[3]]
                        gender_FvKO_SE=modeloutput[["tTable"]][[(3+lengthoftable)]]
                        gender_FvKO_p_value=modeloutput[["tTable"]][[(3+4*lengthoftable)]]
                        gender_MvKO_estimate=modeloutput[["tTable"]][[4]]
                        gender_MvKO_SE=modeloutput[["tTable"]][[(4+lengthoftable)]]
                        gender_MvKO_p_value=modeloutput[["tTable"]][[(4+4*lengthoftable)]]
                        
                    } else if( !Gender_sig && !Interaction_sig){
                        genotype_estimate = modeloutput[["tTable"]][[2]]
                        genotype_estimate_SE = modeloutput[["tTable"]][[(2+lengthoftable)]]
                        genotype_p_value =  modeloutput[["tTable"]][[(2+4*lengthoftable)]]
                        gender_estimate=NA
                        gender_estimate_SE=NA
                        gender_p_value= NA
                        gender_FvKO_estimate= NA
                        gender_FvKO_SE=NA
                        gender_FvKO_p_value=NA
                        gender_MvKO_estimate=NA
                        gender_MvKO_SE=NA
                        gender_MvKO_p_value=NA
                        
                    }else{
                        genotype_estimate = modeloutput[["tTable"]][[2]]
                        genotype_estimate_SE = modeloutput[["tTable"]][[(2+lengthoftable)]]
                        genotype_p_value =  modeloutput[["tTable"]][[(2+4*lengthoftable)]]
                        gender_estimate=modeloutput[["tTable"]][[3]]
                        gender_estimate_SE=modeloutput[["tTable"]][[(3+lengthoftable)]]
                        gender_p_value= modeloutput[["tTable"]][[(3+4*lengthoftable)]]    
                        gender_FvKO_estimate= NA
                        gender_FvKO_SE=NA
                        gender_FvKO_p_value=NA
                        gender_MvKO_estimate=NA
                        gender_MvKO_SE=NA
                        gender_MvKO_p_value=NA
                    }
                    c("Eq1",depVariable, keep_batch, variance_test, Nulltest_genotype_pvalue, genotype_estimate, genotype_estimate_SE,  genotype_p_value,gender_estimate, gender_estimate_SE,  gender_p_value, weight_estimate, weight_estimate_SE, weight_p_value, MM_fitquality, intercept_estimate, intercept_estimate_SE, Interaction_sig, Interaction_test, gender_FvKO_estimate, gender_FvKO_SE,   gender_FvKO_p_value,  gender_MvKO_estimate,  gender_MvKO_SE,  gender_MvKO_p_value)
                    
                }else{
                    #adaption for being a linear model rather than a mixed model
                    intercept_estimate = modeloutput[["tTable"]][[1]]
                    intercept_estimate_SE = modeloutput[["tTable"]][[(1+lengthoftable)]]
                    weight_estimate=NA
                    weight_estimate_SE=NA
                    weight_p_value=NA
                    
                    if((Gender_sig && Interaction_sig) |( !Gender_sig&& Interaction_sig)){
                        
                        genotype_estimate = NA
                        genotype_estimate_SE = NA
                        genotype_p_value =  NA
                        gender_estimate=modeloutput[["tTable"]][[3]]
                        gender_estimate_SE=modeloutput[["tTable"]][[(3+lengthoftable)]]
                        gender_p_value= modeloutput[["tTable"]][[(3+3*lengthoftable)]]
                        gender_FvKO_estimate= modeloutput[["tTable"]][[3]]
                        gender_FvKO_SE=modeloutput[["tTable"]][[(3+lengthoftable)]]
                        gender_FvKO_p_value=modeloutput[["tTable"]][[(3+3*lengthoftable)]]
                        gender_MvKO_estimate=modeloutput[["tTable"]][[4]]
                        gender_MvKO_SE=modeloutput[["tTable"]][[(4+lengthoftable)]]
                        gender_MvKO_p_value=modeloutput[["tTable"]][[(4+3*lengthoftable)]]
                        
                        
                    } else if( !Gender_sig && !Interaction_sig){
                        
                        genotype_estimate = modeloutput[["tTable"]][[2]]
                        genotype_estimate_SE = modeloutput[["tTable"]][[(2+lengthoftable)]]
                        genotype_p_value =  modeloutput[["tTable"]][[(2+3*lengthoftable)]]
                        gender_estimate=NA
                        gender_estimate_SE=NA
                        gender_p_value= NA
                        gender_FvKO_estimate= NA
                        gender_FvKO_SE=NA
                        gender_FvKO_p_value=NA
                        gender_MvKO_estimate=NA
                        gender_MvKO_SE=NA
                        gender_MvKO_p_value=NA
                        
                    }else{
                        
                        genotype_estimate = modeloutput[["tTable"]][[2]]
                        genotype_estimate_SE = modeloutput[["tTable"]][[(2+lengthoftable)]]
                        genotype_p_value =  modeloutput[["tTable"]][[(2+3*lengthoftable)]]
                        gender_estimate=modeloutput[["tTable"]][[3]]
                        gender_estimate_SE=modeloutput[["tTable"]][[(3+lengthoftable)]]
                        gender_p_value= modeloutput[["tTable"]][[(3+3*lengthoftable)]]    
                        gender_FvKO_estimate= NA
                        gender_FvKO_SE=NA
                        gender_FvKO_p_value=NA
                        gender_MvKO_estimate=NA
                        gender_MvKO_SE=NA
                        gender_MvKO_p_value=NA
                    }    
                    c("Eq1",depVariable, keep_batch, variance_test, Nulltest_genotype_pvalue, genotype_estimate, genotype_estimate_SE,  genotype_p_value,gender_estimate, gender_estimate_SE,  gender_p_value, weight_estimate, weight_estimate_SE, weight_p_value, MM_fitquality, intercept_estimate, intercept_estimate_SE, Interaction_sig, Interaction_test, gender_FvKO_estimate, gender_FvKO_SE,   gender_FvKO_p_value,  gender_MvKO_estimate,  gender_MvKO_SE,  gender_MvKO_p_value)
                    
                } 
    },
    withWeight = {
    
                if(!Weight_sig){
                    
                    #If weight is not significant then the output is the same as fitting model Eq1 and so no output is needed. 
                    keep_batch=NA
                    variance_test=NA
                    Nulltest_genotype_pvalue=NA
                    genotype_estimate = NA
                    genotype_estimate_SE = NA
                    genotype_p_value =  NA
                    gender_estimate=NA
                    gender_estimate_SE=NA
                    gender_p_value= NA
                    weight_estimate=NA
                    weight_estimate_SE=NA
                    weight_p_value=NA
                    intercept_estimate = NA
                    intercept_estimate_SE = NA
                    gender_FvKO_estimate= NA
                    gender_FvKO_SE=NA
                    gender_FvKO_p_value=NA
                    gender_MvKO_estimate=NA
                    gender_MvKO_SE=NA
                    gender_MvKO_p_value=NA
                    
                    c("Eq2",depVariable, keep_batch, variance_test, Nulltest_genotype_pvalue, genotype_estimate, genotype_estimate_SE,  genotype_p_value,gender_estimate, gender_estimate_SE,  gender_p_value, weight_estimate, weight_estimate_SE, weight_p_value, MM_fitquality, intercept_estimate, intercept_estimate_SE, Interaction_sig, Interaction_test, gender_FvKO_estimate, gender_FvKO_SE,   gender_FvKO_p_value,  gender_MvKO_estimate,  gender_MvKO_SE,  gender_MvKO_p_value)
                    
                }else{
                    
                    if(keep_batch){
                        
                        #for mixed model 
                        intercept_estimate = modeloutput[["tTable"]][[1]]
                        intercept_estimate_SE = modeloutput[["tTable"]][[(1+lengthoftable)]]
                        
                        if((Weight_sig && Gender_sig && Interaction_sig) | (Weight_sig && !Gender_sig&& Interaction_sig)){
                            
                            genotype_estimate = NA
                            genotype_estimate_SE = NA
                            genotype_p_value =  NA
                            gender_estimate=modeloutput[["tTable"]][[2]]
                            gender_estimate_SE=modeloutput[["tTable"]][[(2+lengthoftable)]]
                            gender_p_value= modeloutput[["tTable"]][[(2+4*lengthoftable)]]
                            gender_FvKO_estimate= modeloutput[["tTable"]][[4]]
                            gender_FvKO_SE=modeloutput[["tTable"]][[(4+lengthoftable)]]
                            gender_FvKO_p_value=modeloutput[["tTable"]][[(4+4*lengthoftable)]]
                            gender_MvKO_estimate=modeloutput[["tTable"]][[5]]
                            gender_MvKO_SE=modeloutput[["tTable"]][[(5+lengthoftable)]]
                            gender_MvKO_p_value=modeloutput[["tTable"]][[(5+4*lengthoftable)]]
                            weight_estimate=modeloutput[["tTable"]][[3]]
                            weight_estimate_SE=modeloutput[["tTable"]][[(3+lengthoftable)]]
                            weight_p_value=modeloutput[["tTable"]][[(3+4*lengthoftable)]]    
                            
                        } else if (Weight_sig && !Gender_sig && !Interaction_sig){    
                            
                            genotype_estimate = modeloutput[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput[["tTable"]][[(2+4*lengthoftable)]]
                            gender_estimate=NA
                            gender_estimate_SE=NA
                            gender_p_value=NA
                            gender_FvKO_estimate= NA
                            gender_FvKO_SE=NA
                            gender_FvKO_p_value=NA
                            gender_MvKO_estimate=NA
                            gender_MvKO_SE=NA
                            gender_MvKO_p_value=NA
                            weight_estimate=modeloutput[["tTable"]][[3]]
                            weight_estimate_SE=modeloutput[["tTable"]][[(3+lengthoftable)]]
                            weight_p_value=modeloutput[["tTable"]][[(3+4*lengthoftable)]]
                            
                        }else if (Weight_sig && Gender_sig && !Interaction_sig){
                            
                            genotype_estimate = modeloutput[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput[["tTable"]][[(2+4*lengthoftable)]]
                            gender_estimate=modeloutput[["tTable"]][[3]]
                            gender_estimate_SE=modeloutput[["tTable"]][[(3+lengthoftable)]]
                            gender_p_value= modeloutput[["tTable"]][[(3+4*lengthoftable)]]    
                            weight_estimate=modeloutput[["tTable"]][[4]]
                            weight_estimate_SE=modeloutput[["tTable"]][[(4+lengthoftable)]]
                            weight_p_value=modeloutput[["tTable"]][[(4+4*lengthoftable)]]
                            gender_FvKO_estimate= NA
                            gender_FvKO_SE=NA
                            gender_FvKO_p_value=NA
                            gender_MvKO_estimate=NA
                            gender_MvKO_SE=NA
                            gender_MvKO_p_value=NA
                            
                        }
                        
                        c("Eq2",depVariable, keep_batch, variance_test, Nulltest_genotype_pvalue, genotype_estimate, genotype_estimate_SE,  genotype_p_value,gender_estimate, gender_estimate_SE,  gender_p_value, weight_estimate, weight_estimate_SE, weight_p_value, MM_fitquality, intercept_estimate, intercept_estimate_SE, Interaction_sig, Interaction_test, gender_FvKO_estimate, gender_FvKO_SE,   gender_FvKO_p_value,  gender_MvKO_estimate,  gender_MvKO_SE,  gender_MvKO_p_value)
                        
                    }else{
                        #adaption for being a linear model rather than a mixed model
                        intercept_estimate = modeloutput[["tTable"]][[1]]
                        intercept_estimate_SE = modeloutput[["tTable"]][[(1+lengthoftable)]]
                        
                        if((Weight_sig && Gender_sig && Interaction_sig )|(Weight_sig && !Gender_sig&& Interaction_sig)){
                            genotype_estimate = NA
                            genotype_estimate_SE = NA
                            genotype_p_value =  NA
                            gender_estimate=modeloutput[["tTable"]][[2]]
                            gender_estimate_SE=modeloutput[["tTable"]][[(2+lengthoftable)]]
                            gender_p_value= modeloutput[["tTable"]][[(2+3*lengthoftable)]]
                            weight_estimate=modeloutput[["tTable"]][[3]]
                            weight_estimate_SE=modeloutput[["tTable"]][[(3+lengthoftable)]]
                            weight_p_value=modeloutput[["tTable"]][[(3+3*lengthoftable)]]    
                            gender_FvKO_estimate= modeloutput[["tTable"]][[4]]
                            gender_FvKO_SE=modeloutput[["tTable"]][[(4+lengthoftable)]]
                            gender_FvKO_p_value=modeloutput[["tTable"]][[(4+3*lengthoftable)]]
                            gender_MvKO_estimate=modeloutput[["tTable"]][[5]]
                            gender_MvKO_SE=modeloutput[["tTable"]][[(5+lengthoftable)]]
                            gender_MvKO_p_value=modeloutput[["tTable"]][[(5+3*lengthoftable)]]        
                            
                        } else if (Weight_sig && Gender_sig && !Interaction_sig){
                            
                            genotype_estimate = modeloutput[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput[["tTable"]][[(2+3*lengthoftable)]]
                            gender_estimate=modeloutput[["tTable"]][[3]]
                            gender_estimate_SE=modeloutput[["tTable"]][[(3+lengthoftable)]]
                            gender_p_value= modeloutput[["tTable"]][[(3+3*lengthoftable)]]
                            weight_estimate=modeloutput[["tTable"]][[4]]
                            weight_estimate_SE=modeloutput[["tTable"]][[(4+lengthoftable)]]
                            weight_p_value=modeloutput[["tTable"]][[(4+3*lengthoftable)]]
                            gender_FvKO_estimate= NA
                            gender_FvKO_SE=NA
                            gender_FvKO_p_value=NA
                            gender_MvKO_estimate=NA
                            gender_MvKO_SE=NA
                            gender_MvKO_p_value=NA
                            
                            
                        }else if (Weight_sig && !Gender_sig && !Interaction_sig){    
                            genotype_estimate = modeloutput[["tTable"]][[2]]
                            genotype_estimate_SE = modeloutput[["tTable"]][[(2+lengthoftable)]]
                            genotype_p_value =  modeloutput[["tTable"]][[(2+3*lengthoftable)]]
                            gender_estimate=NA
                            gender_estimate_SE=NA
                            gender_p_value=NA
                            weight_estimate=modeloutput[["tTable"]][[3]]
                            weight_estimate_SE=modeloutput[["tTable"]][[(3+lengthoftable)]]
                            weight_p_value=modeloutput[["tTable"]][[(3+3*lengthoftable)]]
                            interaction_estimate=NA
                            interaction_estimate_SE=NA
                            interaction_p_value=NA
                            gender_FvKO_estimate= NA
                            gender_FvKO_SE=NA
                            gender_FvKO_p_value=NA
                            gender_MvKO_estimate=NA
                            gender_MvKO_SE=NA
                            gender_MvKO_p_value=NA
                        }    
                        
                        c("Eq2",depVariable, keep_batch, variance_test, Nulltest_genotype_pvalue, genotype_estimate, genotype_estimate_SE,  genotype_p_value,gender_estimate, gender_estimate_SE,  gender_p_value, weight_estimate, weight_estimate_SE, weight_p_value, MM_fitquality, intercept_estimate, intercept_estimate_SE, Interaction_sig, Interaction_test, gender_FvKO_estimate, gender_FvKO_SE,   gender_FvKO_p_value,  gender_MvKO_estimate,  gender_MvKO_SE,  gender_MvKO_p_value)
                    }      
        
                }
       }
      )
    return(values)
}

outputLength <- function(result)
# The output in the vector form depends on effects that are included into the model
{
    if(!is(result,"PhenTestResult")) {
        stop("Please provide result paramtere as PhenTestResult object")
    }
    
    numberofgenders=result$numberGenders
    keep_weight <- result$weightEffect
    keep_gender <- result$genderEffect
    keep_interaction <- result$interactionEffect
    keep_batch <- result$batchEffect
    keep_equalvar <- result$varianceEffect
    equation <- result$equation
    
    table_length <- NA
    
    if (equation=="withWeight"){   
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
    # Eq.1
    else{
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
    
    return (table_length)
}
