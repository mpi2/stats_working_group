## Copyright © 2012-2014 EMBL - European Bioinformatics Institute
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##------------------------------------------------------------------------------
## classification.R contains classificationTag function
##------------------------------------------------------------------------------
## Sexual Dimorphism Classification Tag
classificationTag<-function(phenTestResult, userMode="summaryOutput", 
        phenotypeThreshold=0.01,outputMessages=TRUE)
{
    stop_message <- ""
    ChangeClassification <- ""
    ## check PhenTestResult object
    if(is(phenTestResult,"PhenTestResult")) {
        result<-phenTestResult
        
        analysisResults <- analysisResults(phenTestResult)
        depVariable <- getVariable(phenTestResult)
        
        ## stop function if there are no enough input parameters
        if (method(phenTestResult) %in% c("MM","TF","LR")){
            equation <- analysisResults$equation
            keep_weight <- analysisResults$model.effect.weight
            keep_sex <- analysisResults$model.effect.sex
            keep_interaction <- analysisResults$model.effect.interaction
            keep_batch <- analysisResults$model.effect.batch
            keep_equalvar <- analysisResults$model.effect.variance
            model.output <- analysisResults$model.output
            
            if (is.null(equation) || is.null(depVariable) 
                    || is.null(keep_batch) || is.null(keep_equalvar)
                    || is.null(keep_sex) || is.null(keep_interaction)) {
                stop_message <- "Error:\nPlease run function 'testDataset' first.\n"}
            }
        else
            if (is.null(length(analysisResults) == 0))
                stop_message <- "Error:\nPlease run function 'testDataset' first.\n"
    }
    else{
        stop_message <- "Error:\nPlease create a PhenTestResult object first by using function 'testDataset'.\n"
    }
    
    if (!(userMode %in% c("summaryOutput","vectorOutput"))) {
        stop_message <- paste(stop_message,
                "Error:\nPlease define 'userMode' you would like to use from the following options:",
                " 'summaryOutout' or 'vectorOutput'.\n",sep="")
    }
    
    if (nchar(stop_message)>0){
        if (outputMessages){
            message(stop_message)
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt))
            stop()
        }
        else {
            stop(stop_message)
        }
    }
    else {
        # MM AND TF 
        if (method(phenTestResult) %in% c("MM","TF","LR")) {
            if (userMode=="summaryOutput") {
                if(is.na(analysisResults$model.output.genotype.nulltest.pVal)==TRUE){
                    ChangeClassification <- NA
                }else 
                if(analysisResults$model.output.genotype.nulltest.pVal>phenotypeThreshold){
                    if (analysisResults$numberSexes==1){
                        ChangeClassification <- paste("With phenotype threshold value",
                                phenotypeThreshold,"- no significant change for the one sex tested")
                    }
                    else {
                        ChangeClassification <- paste("With phenotype threshold value",
                                phenotypeThreshold,"- no significant change")
                    }
                }else{
                    if(analysisResults$model.effect.interaction==FALSE) {
                        if (analysisResults$numberSexes==1){
                            ChangeClassification <- paste("With phenotype threshold value",phenotypeThreshold,
                                    "- a significant change for the one sex tested")
                        }
                        else{
                            ChangeClassification <- paste("With phenotype threshold value",
                                    phenotypeThreshold,"- both sexes equally")
                        }
                    }else 
                    if(analysisResults$model.output.summary["sex_FvKO_p_value"]>=0.05 
                            && analysisResults$model.output.summary["sex_MvKO_p_value"]>=0.05){
                        ChangeClassification <- paste("With phenotype threshold value",
                                phenotypeThreshold,"- cannot classify effect")
                    }else 
                    if(analysisResults$model.output.summary["sex_FvKO_p_value"]<0.05 
                            && analysisResults$model.output.summary["sex_MvKO_p_value"]>=0.05){
                        ChangeClassification <- paste("With phenotype threshold value",
                                phenotypeThreshold,"- females only")
                    }else 
                    if(analysisResults$model.output.summary["sex_FvKO_p_value"]>=0.05 
                            && analysisResults$model.output.summary["sex_MvKO_p_value"]<0.05){
                        ChangeClassification <- paste("With phenotype threshold value",
                                phenotypeThreshold,"- males only")
                    }else 
                    if(analysisResults$model.output.summary["sex_FvKO_estimate"]>0 && 
                            analysisResults$model.output.summary["sex_MvKO_estimate"]>0 | 
                            analysisResults$model.output.summary["sex_FvKO_estimate"]<0 && 
                            analysisResults$model.output.summary["sex_MvKO_estimate"]<0){
                        if(abs(analysisResults$model.output.summary["sex_FvKO_estimate"])
                                >abs(analysisResults$model.output.summary["sex_MvKO_estimate"])){
                            ChangeClassification <- paste("With phenotype threshold value",phenotypeThreshold,
                                    "- different size as females greater")
                            ## change could be positive or negative but size change greater
                        }else{
                            ChangeClassification <- paste("With phenotype threshold value",
                                    phenotypeThreshold,"- different size as males greater")
                        }
                    }else{
                        ChangeClassification <- paste("With phenotype threshold value",
                                phenotypeThreshold,"- different direction for the sexes")
                    }
                }
            }
            #vectorOutput
            else {
                
                if(analysisResults$model.effect.interaction==FALSE) {
                    if (analysisResults$numberSexes==1){
                        ChangeClassification <- paste("If phenotype is significant it is for the one sex tested")
                    }else{
                        ChangeClassification <- paste("If phenotype is significant - both sexes equally")
                    }
                }else 
                if(analysisResults$model.output.summary["sex_FvKO_p_value"]>=0.05 
                        && analysisResults$model.output.summary["sex_MvKO_p_value"]>=0.05){
                    ChangeClassification <- paste("If phenotype is significant - can not classify effect")
                }else 
                if(analysisResults$model.output.summary["sex_FvKO_p_value"]<0.05 
                        && analysisResults$model.output.summary["sex_MvKO_p_value"]>=0.05){
                    ChangeClassification <- paste("If phenotype is significant - females only")
                }else 
                if(analysisResults$model.output.summary["sex_FvKO_p_value"]>=0.05 
                        && analysisResults$model.output.summary["sex_MvKO_p_value"]<0.05){
                    ChangeClassification <- paste("If phenotype is significant - males only")
                }else 
                if(analysisResults$model.output.summary["sex_FvKO_estimate"]>0 && 
                        analysisResults$model.output.summary["sex_MvKO_estimate"]>0 |
                        analysisResults$model.output.summary["sex_FvKO_estimate"]<0 && 
                        analysisResults$model.output.summary["sex_MvKO_estimate"]<0){
                    if(abs(analysisResults$model.output.summary["sex_FvKO_estimate"])
                            >abs(analysisResults$model.output.summary["sex_MvKO_estimate"])){
                        ChangeClassification <- paste("If phenotype is significant - different size as females greater")
                        ## change could be positive or negative but size change greater
                    }else{
                        ChangeClassification <- paste("If phenotype is significant - different size as males greater")
                    }
                }else{
                    ChangeClassification <- paste("If phenotype is significant - different direction for the sexes")
                }
                
            }
        }
        # FISHER EXACT TEST
        else if (method(phenTestResult) =="FE"){
            if (userMode=="summaryOutput") { 
                male_p.value <- 10
                female_p.value <- 10
                all_p.value <- 10
                for (i in seq_along(analysisResults)) {
                    val <- analysisResults[[i]]
                    if (analysedSubset(val)=="all"){
                        all_p.value<-as.numeric(as.vector(getColumnView(val))[2]) 
                               
                    }
                    if (analysedSubset(val)=="males"){
                        male_p.value<-as.numeric(as.vector(getColumnView(val))[2])             
                    }
                    if (analysedSubset(val)=="females"){
                        female_p.value<-as.numeric(as.vector(getColumnView(val))[2])                
                    }
                }
                
                if (noSexes(phenTestResult)==1){
                    ChangeClassification <- paste("Not significant for the one sex tested")
                }
                else {
                    ChangeClassification <- paste("Not significant")
                }
                
                
                # Tag
                # combined & males & females
                if(all_p.value < phenotypeThreshold && 
                        male_p.value < phenotypeThreshold && 
                        female_p.value < phenotypeThreshold)
                    ChangeClassification <- paste("With phenotype threshold value",phenotypeThreshold, 
                        "- significant in males, females and in combined dataset")
                # combined & males & !females
                if(all_p.value < phenotypeThreshold && 
                        male_p.value < phenotypeThreshold && 
                        female_p.value >= phenotypeThreshold)
                    ChangeClassification <- paste("With phenotype threshold value",phenotypeThreshold, 
                        "- significant in males and in combined dataset")
                # combined & !males & females
                if(all_p.value < phenotypeThreshold && 
                        male_p.value >= phenotypeThreshold && 
                        female_p.value < phenotypeThreshold)
                    ChangeClassification <- paste("With phenotype threshold value",phenotypeThreshold, 
                        "- significant in females and in combined dataset")
                # combined & !males & !females
                if(all_p.value < phenotypeThreshold && 
                        male_p.value >= phenotypeThreshold && 
                        female_p.value >= phenotypeThreshold){
                    if (noSexes(phenTestResult)==2)
                        ChangeClassification <- paste("With phenotype threshold value",phenotypeThreshold, 
                            "- significant in combined dataset only")
                    else
                        ChangeClassification <- paste("With phenotype threshold value",phenotypeThreshold, 
                            "- significant for the sex tested")
                }
                # !combined & males & females
                if(all_p.value >= phenotypeThreshold && 
                        male_p.value < phenotypeThreshold && 
                        female_p.value < phenotypeThreshold)
                    ChangeClassification <- paste("With phenotype threshold value",phenotypeThreshold, 
                        "- significant in males and in females datasets")
                # !combined & males & !females
                if(all_p.value >= phenotypeThreshold && 
                        male_p.value < phenotypeThreshold && 
                        female_p.value >= phenotypeThreshold)
                    ChangeClassification <- paste("With phenotype threshold value",phenotypeThreshold, 
                        "- significant in males dataset only")
                # !combined & !males & females
                if(all_p.value >= phenotypeThreshold && 
                        male_p.value >= phenotypeThreshold && 
                        female_p.value < phenotypeThreshold)
                    ChangeClassification <- paste("With phenotype threshold value",phenotypeThreshold, 
                        "- significant in females dataset only")
          }
          else {
                ChangeClassification <- "NA"
          }
        }
        # REFERENCE RANGE PLUS
        else if (method(phenTestResult)=="RR"){
            if (userMode=="summaryOutput") { 
                direction_all <- "NA"
                all_p.value <- 10
                direction_females <- "NA"
                female_p.value <- 10
                direction_males <- "NA"
                male_p.value <- 10
                high_male_p.value <- 10
                high_female_p.value <- 10
                high_all_p.value <- 10
                low_male_p.value <- 10
                low_female_p.value <- 10
                low_all_p.value <- 10
                
                for (i in seq_along(analysisResults)) {
                    val <- analysisResults[[i]]
                    if (comparison(val)=="High vs Normal/Low") {
                        if (analysedSubset(val)=="all"){
                            high_all_p.value<-as.numeric(as.vector(getColumnView(val))[2])      
                        }
                        if (analysedSubset(val)=="males"){
                            high_male_p.value<-as.numeric(as.vector(getColumnView(val))[2])                  
                        }
                        if (analysedSubset(val)=="females"){
                            high_female_p.value<-as.numeric(as.vector(getColumnView(val))[2])                  
                        }
                    }
                    else {
                        if (analysedSubset(val)=="all"){
                            low_all_p.value<-as.numeric(as.vector(getColumnView(val))[2])                  
                        }
                        if (analysedSubset(val)=="males"){
                            low_male_p.value<-as.numeric(as.vector(getColumnView(val))[2])                  
                        }
                        if (analysedSubset(val)=="females"){
                            low_female_p.value<-as.numeric(as.vector(getColumnView(val))[2])                  
                        }
                    }
                }
                
                # High classification p-val is less than threshold and low classification p-val is more that threshold
                if (high_all_p.value < phenotypeThreshold && low_all_p.value >= phenotypeThreshold){
                    direction_all <- "High"
                    all_p.value <- high_all_p.value
                }
                # Low classification p-val is less than threshold and high classification p-val is more that threshold
                else if (high_all_p.value >= phenotypeThreshold && low_all_p.value < phenotypeThreshold){
                    direction_all <- "Low"
                    all_p.value <- low_all_p.value
                }
                
                    if (high_male_p.value < phenotypeThreshold && low_male_p.value >= phenotypeThreshold){
                        direction_males <- "High"
                        male_p.value <- high_male_p.value
                    }
                    else if (high_male_p.value >= phenotypeThreshold && low_male_p.value < phenotypeThreshold){
                        direction_males <- "Low"
                        male_p.value <- low_male_p.value
                    }
   
                
                    if (high_female_p.value < phenotypeThreshold && low_female_p.value >= phenotypeThreshold){
                        direction_females <- "High"
                        female_p.value <- high_female_p.value
                    }
                    else if (high_female_p.value >= phenotypeThreshold && low_female_p.value < phenotypeThreshold){
                        direction_females <- "Low"
                        female_p.value <- low_female_p.value
                    }
 
                
                if (noSexes(phenTestResult)==1){
                    ChangeClassification <- paste("Not significant for the one sex tested")
                }
                else {
                    ChangeClassification <- paste("Not significant")
                }
                    
                if (noSexes(phenTestResult)==2) {
                    # Tag
                    # combined & males & females
                    if(all_p.value < phenotypeThreshold && male_p.value < phenotypeThreshold 
                            && female_p.value < phenotypeThreshold)
                    ChangeClassification <- paste("With phenotype threshold value ",phenotypeThreshold,
                                    " - significant in males (",direction_males,
                                    "), females (",direction_females,") and in combined dataset (",direction_all,")",sep="")
                    # combined & males & !females
                    if(all_p.value < phenotypeThreshold && male_p.value < phenotypeThreshold 
                            && female_p.value >= phenotypeThreshold)
                    ChangeClassification <- paste("With phenotype threshold value ",phenotypeThreshold,
                                    " - significant in males (",direction_males,
                                    ") and in combined dataset (",direction_all,")",sep="")
                    # combined & !males & females
                    if(all_p.value < phenotypeThreshold && male_p.value >= phenotypeThreshold 
                            && female_p.value < phenotypeThreshold)
                    ChangeClassification <- paste("With phenotype threshold value ",phenotypeThreshold,
                                    " - significant in females (",direction_females,
                                    ") and in combined dataset (",direction_all,")",sep="")
                    # combined & !males & !females
                    if(all_p.value < phenotypeThreshold && male_p.value >= phenotypeThreshold 
                            && female_p.value >= phenotypeThreshold){
                            ChangeClassification <- paste("With phenotype threshold value ",phenotypeThreshold, 
                                " - significant in combined dataset only (",direction_all,")",sep="")
                            
                    }
                    # !combined & males & females
                    if(all_p.value >= phenotypeThreshold && male_p.value < phenotypeThreshold 
                            && female_p.value < phenotypeThreshold)
                        ChangeClassification <- paste("With phenotype threshold value ",phenotypeThreshold,
                                    " - significant in males (",direction_males,
                                    ") and females (",direction_females,") datasets",sep="") 
                    # !combined & males & !females
                    if(all_p.value >= phenotypeThreshold && male_p.value < phenotypeThreshold
                            && female_p.value >= phenotypeThreshold)
                        ChangeClassification <- paste("With phenotype threshold value ",phenotypeThreshold,
                            " - significant in males (",direction_males,") dataset only",sep="")
                    # !combined & !males & females
                    if(all_p.value >= phenotypeThreshold && male_p.value >= phenotypeThreshold 
                            && female_p.value < phenotypeThreshold)
                        ChangeClassification <- paste("With phenotype threshold value ",phenotypeThreshold,
                            " - significant in females (",direction_females,") dataset only",sep="")
                }
                else {
                    if(all_p.value < phenotypeThreshold) {
                        ChangeClassification <- paste("With phenotype threshold value ",phenotypeThreshold,
                                " - significant for the sex tested (",direction_all,")",sep="")
                    }
                }
            }
            else {
                ChangeClassification <- "NA"
            }
        }
        return(ChangeClassification)
        
    }
}
##------------------------------------------------------------------------------