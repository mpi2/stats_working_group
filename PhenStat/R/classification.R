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
    genderVector<-vectorOutput(result)[25:30] 
    m<-matrix(genderVector,1,6)
    colnames(m) = c("FvKO_par", "FvKO_SE", "FvKO_pvalue", 
            "MvKO_par", "MvKO_SE", "MvKO_pvalue")
    gvr<-as.data.frame(m,stringsAsFactors=FALSE)
    
        
        #added a layer for no fitting of data ie NA
        if(is.na(result$model.output.genotype.nulltest.pVal)==TRUE){
            ChangeClassification==NA
        }else if(result$model.output.genotype.nulltest.pVal>phenotypeThreshold){
            ChangeClassification=paste("With phenotype threshold value",phenotypeThreshold,"- no significant change")
        }else{
            if(result$model.effect.interaction==FALSE) {
                ChangeClassification=paste("With phenotype threshold value",phenotypeThreshold,"- both sexes equally")
            
            }else if(gvr$FvKO_pvalue>=0.05 && gvrMvKO_pvalue>=0.05){
                ChangeClassification=paste("With phenotype threshold value",phenotypeThreshold,"- cannot classify effect")
            }else if(gvr$FvKO_pvalue<0.05 && gvr$MvKO_pvalue>=0.05){
                ChangeClassification=paste("With phenotype threshold value",phenotypeThreshold,"- females only")
            }else if(gvr$FvKO_pvalue>=0.05 && gvr$MvKO_pvalue<0.05){
                ChangeClassificationpaste("With phenotype threshold value",phenotypeThreshold,"- males only")
            }else if(gvr$FvKO_par>0 && gvr$MvKO_par>0 |gvr$FvKO_par<0 && gvr$MvKO_par<0){
                if(abs(gvr$FvKO_par)>abs(gvr$MvKO_par)){
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

