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
#-----------------------------------------------------------------------------------
# vectorOutput.R contains vectorOutput and outputLength functions
#-----------------------------------------------------------------------------------
vectorOutput <- function(phenTestResult)
# Wrapper to prepare the output of the modeling and testing results in vector form. Assumes that modeling results are 
# stored in the phenTestResult object (output from functions testDataset and buildFinalModel)
{
    
    if (phenTestResult$method=="MM") {
        equation <- switch(phenTestResult$equation,withoutWeight = {"Eq1"},withWeight = {"Eq2"})
        
        classificationValue <- classificationTag(phenTestResult,userMode="vectorOutput")
        
        vectorOutput <- c(equation,phenTestResult$depVariable, phenTestResult$model.effect.batch, 
                phenTestResult$model.effect.variance, 
                phenTestResult$model.output.genotype.nulltest.pVal, 
                phenTestResult$model.output.summary["genotype_estimate"], 
                phenTestResult$model.output.summary["genotype_estimate_SE"],  
                phenTestResult$model.output.summary["genotype_p_value"],
                phenTestResult$model.output.summary["gender_estimate"], 
                phenTestResult$model.output.summary["gender_estimate_SE"],  
                phenTestResult$model.output.summary["gender_p_value"], 
                phenTestResult$model.output.summary["weight_estimate"], 
                phenTestResult$model.output.summary["weight_estimate_SE"], 
                phenTestResult$model.output.summary["weight_p_value"], 
                phenTestResult$model.output.quality, 
                phenTestResult$model.output.summary["intercept_estimate"], 
                phenTestResult$model.output.summary["intercept_estimate_SE"], 
                phenTestResult$model.output.interaction,
                phenTestResult$model.output.summary["gender_FvKO_estimate"], 
                phenTestResult$model.output.summary["gender_FvKO_SE"], 
                phenTestResult$model.output.summary["gender_FvKO_p_value"],  
                phenTestResult$model.output.summary["gender_MvKO_estimate"],  
                phenTestResult$model.output.summary["gender_MvKO_SE"], 
                phenTestResult$model.output.summary["gender_MvKO_p_value"],
                classificationValue)
    }
    else if (phenTestResult$method=="FE"){
        male_pval <- NA
        female_pval <- NA
        if (!is.null(phenTestResult$model.output$male)){
            male_pval<-round(phenTestResult$model.output$male$p.val,digits=3)
        }
        if (!is.null(phenTestResult$model.output$female)){
            female_pval<-round(phenTestResult$model.output$female$p.val,digits=3)
        }
        
        vectorOutput <- c("Fisher Exact Test",phenTestResult$depVariable, phenTestResult$model.effect.batch, 
                phenTestResult$model.effect.variance,
                round(phenTestResult$model.output$all$p.val,digits=3),
                NA, 
                NA,  
                NA,
                NA, 
                NA,  
                NA, 
                NA, 
                NA, 
                NA, 
                NA, 
                NA, 
                NA, 
                NA,
                NA, 
                NA, 
                female_pval,  
                NA,  
                NA, 
                male_pval,
                NA)
    }
    names(vectorOutput)<-NULL
    return(vectorOutput)
}