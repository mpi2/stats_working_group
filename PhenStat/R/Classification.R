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

# classification.R contains classificationTag function

classificationTag<-function(result, phenotypeThreshold=0.01)
# Sexual Dimorphism Classification Tag

{
        #added a layer for no fitting of data ie NA
        if(is.na(result$model.output.genotype.nulltest.pVal)==TRUE){
            ChangeClassification==NA
        }else if(result$model.output.genotype.nulltest.pVal>phenotypeThreshold){
            ChangeClassification=paste("With phenotype threshold value",phenotypeThreshold,"- no significant change")
        }else{
            if(result$model.effect.interaction==FALSE) {
                ChangeClassification=paste("With phenotype threshold value",phenotypeThreshold,"- both sexes equally")
            
            }else if(result$model.output.summary["gender_FvKO_p_value"]>=0.05 && result$model.output.summary["gender_MvKO_p_value"]>=0.05){
                ChangeClassification=paste("With phenotype threshold value",phenotypeThreshold,"- cannot classify effect")
            }else if(result$model.output.summary["gender_FvKO_p_value"]<0.05 && result$model.output.summary["gender_MvKO_p_value"]>=0.05){
                ChangeClassification=paste("With phenotype threshold value",phenotypeThreshold,"- females only")
            }else if(result$model.output.summary["gender_FvKO_p_value"]>=0.05 && result$model.output.summary["gender_MvKO_p_value"]<0.05){
                ChangeClassificationpaste("With phenotype threshold value",phenotypeThreshold,"- males only")
            }else if(result$model.output.summary["gender_FvKO_estimate"]>0 && result$model.output.summary["gender_MvKO_estimate"]>0 |result$model.output.summary["gender_FvKO_estimate"]<0 && result$model.output.summary["gender_MvKO_estimate"]<0){
                if(abs(result$model.output.summary["gender_FvKO_estimate"])>abs(result$model.output.summary["gender_MvKO_estimate"])){
                    ChangeClassification=paste("With phenotype threshold value",phenotypeThreshold,"- different size as females greater") 
                    # change could be positive or negative but size change greater
                }else{
                    ChangeClassification=paste("With phenotype threshold value",phenotypeThreshold,"- different size as males greater")
                }        
            }else{
                ChangeClassification=paste("With phenotype threshold value",phenotypeThreshold,"- different direction for the sexes")
            }                                        
        }
    
    return(ChangeClassification)
}

