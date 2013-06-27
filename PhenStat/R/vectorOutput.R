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

# vectorOutput.R contains vectorOutput and outputLength functions

vectorOutput <- function(result)
# Wrapper to prepare the output of the modeling and testing results in vector form
{
    if(!is(result,"PhenTestResult")) {
        stop ("Please provide PhenTestResult object as an argument")
    }
    
    if (is.null(result$model.formula.genotype)) stop("There are no results to wrap. Please run function 'buildFinalModel' first")
    
    equation <- switch(result$equation,withoutWeight = {"Eq1"},withWeight = {"Eq2"})
    
    vectorOutput <- c(equation,result$depVariable, result$model.effect.batch, result$model.effect.variance, result$model.output.genotype.nulltest.pVal, 
            result$model.output.summary["genotype_estimate"], 
            result$model.output.summary["genotype_estimate_SE"],  
            result$model.output.summary["genotype_p_value"],
            result$model.output.summary["gender_estimate"], 
            result$model.output.summary["gender_estimate_SE"],  
            result$model.output.summary["gender_p_value"], 
            result$model.output.summary["weight_estimate"], 
            result$model.output.summary["weight_estimate_SE"], 
            result$model.output.summary["weight_p_value"], 
            result$model.output.quality, 
            result$model.output.summary["intercept_estimate"], 
            result$model.output.summary["intercept_estimate_SE"], 
            result$model.output.interaction,
            result$model.output.summary["gender_FvKO_estimate"], 
            result$model.output.summary["gender_FvKO_SE"], 
            result$model.output.summary["gender_FvKO_p_value"],  
            result$model.output.summary["gender_MvKO_estimate"],  
            result$model.output.summary["gender_MvKO_SE"], 
            result$model.output.summary["gender_MvKO_p_value"])
    names(vectorOutput)<-NULL
    return(vectorOutput)
}

outputLength <- function(result)
# The output in the vector form depends on effects that are included into the model
{
    if(!is(result,"PhenTestResult")) {
        stop("Please provide result paramtere as PhenTestResult object")
    }
    
    numberofgenders=result$numberGenders
    keep_weight <- result$model.effect.weight
    keep_gender <- result$model.effect.gender
    keep_interaction <- result$model.effect.interaction
    keep_batch <- result$model.effect.batch
    keep_equalvar <- result$model.effect.variance
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
