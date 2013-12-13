## Copyright © 2011-2013 EMBL - European Bioinformatics Institute
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
    ## check PhenTestResult object
    if(is(phenTestResult,"PhenTestResult")) {
        result<-phenTestResult
        depVariable <- result$depVariable
        equation <- result$equation
        keep_weight <- result$model.effect.weight
        keep_gender <- result$model.effect.gender
        keep_interaction <- result$model.effect.interaction
        keep_batch <- result$model.effect.batch
        keep_equalvar <- result$model.effect.variance
        model.output <- result$model.output
        
        ## stop function if there are no enough input parameters
        if (phenTestResult$method=="MM"){
            if (is.null(equation) || is.null(depVariable) 
                    || is.null(keep_batch) || is.null(keep_equalvar)
                    || is.null(keep_gender) || is.null(keep_interaction)) {
                stop_message <- "Error:\nPlease run function 'testDataset' first
                .\n"}}
        else
        if (is.null(model.output))
        stop_message <- "Error:\nPlease run function 'testDataset' first.\n"
    }
    else{
        stop_message <- "Error:\nPlease create a PhenTestResult object first by 
        using function 'testDataset'.\n"
    }
    
    if (!(userMode %in% c("summaryOutput","vectorOutput"))) {
        stop_message <- paste(stop_message,"Error:\nPlease define 'userMode' you 
                would like to use from the following options: 'summaryOutout' 
                or 'vectorOutput'.\n",sep="")
    }
    
    if (nchar(stop_message)>0){
        if (outputMessages){
            message(stop_message)
        }
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))
        stop()
    }
    else {
        if (phenTestResult$method=="MM") {
            if (userMode=="summaryOutput") {
                if(is.na(result$model.output.genotype.nulltest.pVal)==TRUE){
                    ChangeClassification==NA
                }else 
                if(result$model.output.genotype.nulltest.pVal>phenotypeThreshold){
                    ChangeClassification <- paste("With phenotype threshold 
                            value",phenotypeThreshold,"- no significant change")
                }else{
                    if(result$model.effect.interaction==FALSE) {
                        if (result$numberGenders==1){
                            ChangeClassification <- paste("With phenotype 
                                    threshold value",phenotypeThreshold,
                                    "- a significant change for the one genotype tested")
                        }
                        else{
                            ChangeClassification <- paste("With phenotype 
                                    threshold value",phenotypeThreshold,"- both sexes equally")
                        }
                    }else 
                    if(result$model.output.summary["gender_FvKO_p_value"]>=0.05 
                            && result$model.output.summary["gender_MvKO_p_value"]>=0.05){
                        ChangeClassification <- paste("With phenotype threshold 
                                value",phenotypeThreshold,"- cannot classify effect")
                    }else 
                    if(result$model.output.summary["gender_FvKO_p_value"]<0.05 
                            && result$model.output.summary["gender_MvKO_p_value"]>=0.05){
                        ChangeClassification <- paste("With phenotype threshold 
                                value",phenotypeThreshold,"- females only")
                    }else 
                    if(result$model.output.summary["gender_FvKO_p_value"]>=0.05 
                            && result$model.output.summary["gender_MvKO_p_value"]<0.05){
                        ChangeClassification <- paste("With phenotype threshold 
                                value",phenotypeThreshold,"- males only")
                    }else 
                    if(result$model.output.summary["gender_FvKO_estimate"]>0 && 
                            result$model.output.summary["gender_MvKO_estimate"]>0 | 
                            result$model.output.summary["gender_FvKO_estimate"]<0 && 
                            result$model.output.summary["gender_MvKO_estimate"]<0){
                        if(abs(result$model.output.summary["gender_FvKO_estimate"])
                                >abs(result$model.output.summary["gender_MvKO_estimate"])){
                            ChangeClassification <- paste("With phenotype threshold 
                                    value",phenotypeThreshold,
                                    "- different size as females greater")
                            ## change could be positive or negative but size change greater
                        }else{
                            ChangeClassification <- paste("With phenotype 
                                    threshold value",phenotypeThreshold,
                                    "- different size as males greater")
                        }
                    }else{
                        ChangeClassification <- paste("With phenotype threshold 
                                value",phenotypeThreshold,
                                "- different direction for the sexes")
                    }
                }
            }
            else {
                
                if(result$model.effect.interaction==FALSE) {
                    if (result$numberGenders==1){
                        ChangeClassification <- paste("If phenotype is significant 
                                it is for the one genotype tested")
                    }else{
                        ChangeClassification <- paste("If phenotype is significant 
                                - both sexes equally")
                    }
                }else 
                if(result$model.output.summary["gender_FvKO_p_value"]>=0.05 
                        && result$model.output.summary["gender_MvKO_p_value"]>=0.05){
                    ChangeClassification <- paste("If phenotype is significant - 
                            cannot classify effect")
                }else 
                if(result$model.output.summary["gender_FvKO_p_value"]<0.05 
                        && result$model.output.summary["gender_MvKO_p_value"]>=0.05){
                    ChangeClassification <- paste("If phenotype is significant - 
                            females only")
                }else 
                if(result$model.output.summary["gender_FvKO_p_value"]>=0.05 
                        && result$model.output.summary["gender_MvKO_p_value"]<0.05){
                    ChangeClassification <- paste("If phenotype is significant - 
                            males only")
                }else 
                if(result$model.output.summary["gender_FvKO_estimate"]>0 && 
                        result$model.output.summary["gender_MvKO_estimate"]>0 |
                        result$model.output.summary["gender_FvKO_estimate"]<0 && 
                        result$model.output.summary["gender_MvKO_estimate"]<0){
                    if(abs(result$model.output.summary["gender_FvKO_estimate"])
                            >abs(result$model.output.summary["gender_MvKO_estimate"])){
                        ChangeClassification <- paste("If phenotype is significant -
                                different size as females greater")
                        ## change could be positive or negative but size change greater
                    }else{
                        ChangeClassification <- paste("If phenotype is significant -
                                different size as males greater")
                    }
                }else{
                    ChangeClassification <- paste("If phenotype is significant - 
                            different direction for the sexes")
                }
                
            }
        }
        else if (phenTestResult$method=="FE"){
            if (!is.null(phenTestResult$model.output$male)){
                male_p.value <- result$model.output$male$p.value
            }
            else {
                male_p.value <- 10
            }
            if (!is.null(phenTestResult$model.output$female)){
                female_p.value <- result$model.output$female$p.value
            }
            else {
                female_p.value <- 10
            }
            all_p.value <- result$model.output$all$p.value
            
            ChangeClassification <- paste("Not significant")
            # Tag
            if(all_p.value < 0.05 && male_p.value < 0.05 
                    && female_p.value < 0.05)
            ChangeClassification <- paste("Significant in males, 
                    females and in combined dataset")
            
            if(all_p.value < 0.05 && male_p.value < 0.05 
                    && female_p.value >= 0.05)
            ChangeClassification <- paste("Significant in males and in 
                    combined dataset")
            if(all_p.value < 0.05 && male_p.value >= 0.05 
                    && female_p.value >= 0.05)
            ChangeClassification <- paste("Significant in females and in 
                    combined dataset")
            if(all_p.value < 0.05 && male_p.value >= 0.05 
                    && female_p.value >= 0.05){
                if (phenTestResult$numberGenders==2)
                ChangeClassification <- paste("Significant in combined 
                        dataset only")
                else
                ChangeClassification <- paste("Significant for the sex 
                        tested")
            }
            if(all_p.value >= 0.05 && male_p.value < 0.05 
                    && female_p.value < 0.05)
            ChangeClassification <- paste("Significant in males and in 
                    females datasets")
            if(all_p.value >= 0.05 && male_p.value < 0.05 
                    && female_p.value >= 0.05)
            ChangeClassification <- paste("Significant in males dataset only")
            if(all_p.value >= 0.05 && male_p.value >= 0.05 
                    && female_p.value < 0.05)
            ChangeClassification <- paste("Significant in females dataset 
                    only")
        }
        return(ChangeClassification)
        
    }
}
##------------------------------------------------------------------------------