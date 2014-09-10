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
        depVariable <- result$depVariable
        equation <- result$equation
        keep_weight <- result$model.effect.weight
        keep_sex <- result$model.effect.sex
        keep_interaction <- result$model.effect.interaction
        keep_batch <- result$model.effect.batch
        keep_equalvar <- result$model.effect.variance
        model.output <- result$model.output
        
        ## stop function if there are no enough input parameters
        if (phenTestResult$method=="MM"){
            if (is.null(equation) || is.null(depVariable) 
                    || is.null(keep_batch) || is.null(keep_equalvar)
                    || is.null(keep_sex) || is.null(keep_interaction)) {
                stop_message <- "Error:\nPlease run function 'testDataset' first.\n"}}
        else
        if (is.null(model.output))
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
        if (phenTestResult$method %in% c("MM","TF")) {
            if (userMode=="summaryOutput") {
                if(is.na(result$model.output.genotype.nulltest.pVal)==TRUE){
                    ChangeClassification <- NA
                }else 
                if(result$model.output.genotype.nulltest.pVal>phenotypeThreshold){
                    ChangeClassification <- paste("With phenotype threshold value",
                            phenotypeThreshold,"- no significant change")
                }else{
                    if(result$model.effect.interaction==FALSE) {
                        if (result$numberSexes==1){
                            ChangeClassification <- paste("With phenotype threshold value",phenotypeThreshold,
                                    "- a significant change for the one sex tested")
                        }
                        else{
                            ChangeClassification <- paste("With phenotype threshold value",
                                    phenotypeThreshold,"- both sexes equally")
                        }
                    }else 
                    if(result$model.output.summary["sex_FvKO_p_value"]>=0.05 
                            && result$model.output.summary["sex_MvKO_p_value"]>=0.05){
                        ChangeClassification <- paste("With phenotype threshold value",
                                phenotypeThreshold,"- cannot classify effect")
                    }else 
                    if(result$model.output.summary["sex_FvKO_p_value"]<0.05 
                            && result$model.output.summary["sex_MvKO_p_value"]>=0.05){
                        ChangeClassification <- paste("With phenotype threshold value",
                                phenotypeThreshold,"- females only")
                    }else 
                    if(result$model.output.summary["sex_FvKO_p_value"]>=0.05 
                            && result$model.output.summary["sex_MvKO_p_value"]<0.05){
                        ChangeClassification <- paste("With phenotype threshold value",
                                phenotypeThreshold,"- males only")
                    }else 
                    if(result$model.output.summary["sex_FvKO_estimate"]>0 && 
                            result$model.output.summary["sex_MvKO_estimate"]>0 | 
                            result$model.output.summary["sex_FvKO_estimate"]<0 && 
                            result$model.output.summary["sex_MvKO_estimate"]<0){
                        if(abs(result$model.output.summary["sex_FvKO_estimate"])
                                >abs(result$model.output.summary["sex_MvKO_estimate"])){
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
                
                if(result$model.effect.interaction==FALSE) {
                    if (result$numberSexes==1){
                        ChangeClassification <- paste("If phenotype is significant it is for the one sex tested")
                    }else{
                        ChangeClassification <- paste("If phenotype is significant - both sexes equally")
                    }
                }else 
                if(result$model.output.summary["sex_FvKO_p_value"]>=0.05 
                        && result$model.output.summary["sex_MvKO_p_value"]>=0.05){
                    ChangeClassification <- paste("If phenotype is significant - can not classify effect")
                }else 
                if(result$model.output.summary["sex_FvKO_p_value"]<0.05 
                        && result$model.output.summary["sex_MvKO_p_value"]>=0.05){
                    ChangeClassification <- paste("If phenotype is significant - females only")
                }else 
                if(result$model.output.summary["sex_FvKO_p_value"]>=0.05 
                        && result$model.output.summary["sex_MvKO_p_value"]<0.05){
                    ChangeClassification <- paste("If phenotype is significant - males only")
                }else 
                if(result$model.output.summary["sex_FvKO_estimate"]>0 && 
                        result$model.output.summary["sex_MvKO_estimate"]>0 |
                        result$model.output.summary["sex_FvKO_estimate"]<0 && 
                        result$model.output.summary["sex_MvKO_estimate"]<0){
                    if(abs(result$model.output.summary["sex_FvKO_estimate"])
                            >abs(result$model.output.summary["sex_MvKO_estimate"])){
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
        else if (phenTestResult$method =="FE"){
            if (userMode=="summaryOutput") { 
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
                    if (phenTestResult$numberSexes==2)
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
        else if (phenTestResult$method=="RR"){
            if (userMode=="summaryOutput") { 
                direction_all <- "NA"
                all_p.value <- 10
                direction_females <- "NA"
                female_p.value <- 10
                direction_males <- "NA"
                male_p.value <- 10
                
                RROutput <- phenTestResult$model.output$all
                # High classification p-val is less than threshold and low classification p-val is more that threshold
                if (RROutput[1] < phenotypeThreshold && RROutput[3] >= phenotypeThreshold){
                    direction_all <- "Low"
                    all_p.value <- RROutput[1]
                }
                # Low classification p-val is less than threshold and high classification p-val is more that threshold
                else if (RROutput[1] >= phenotypeThreshold && RROutput[3] < phenotypeThreshold){
                    direction_all <- "High"
                    all_p.value <- RROutput[3]
                }
                
                #direction_females <- names(which.max(result$model.output$percentage_matrix_female[c(1,3),3]))
                #direction_males <- names(which.max(result$model.output$percentage_matrix_male[c(1,3),3]))
                #direction_all <- names(which.max(result$model.output$percentage_matrix_all[c(1,3),3]))
                
                # Low and High have the same Effect sizes
                #if (phenTestResult$numberSexes==2){
                    #   low_es <- result$model.output$percentage_matrix_female[1,3]
                    #   high_es <- result$model.output$percentage_matrix_female[3,3]
                    #   if (low_es == high_es)
                    #   direction_females <- "NA"
                    #   low_es <- result$model.output$percentage_matrix_male[1,3]
                    #   high_es <- result$model.output$percentage_matrix_male[3,3]
                    #   if (low_es == high_es)
                    #   direction_males <- "NA"
                    #}
                #else {
                    #   low_es <- result$model.output$percentage_matrix_all[1,3]
                    #   high_es <- result$model.output$percentage_matrix_all[3,3]
                    #   if (low_es == high_es)
                    #   direction_all <- "NA"
                    #}
                
                if (!is.null(phenTestResult$model.output$male)){
                    RROutput <- phenTestResult$model.output$male
                    if (RROutput[1] < phenotypeThreshold && RROutput[3] >= phenotypeThreshold){
                        direction_males <- "Low"
                        male_p.value <- RROutput[1]
                    }
                    else if (RROutput[1] >= phenotypeThreshold && RROutput[3] < phenotypeThreshold){
                        direction_males <- "High"
                        male_p.value <- RROutput[3]
                    }
                    #male_p.value <- result$model.output$male$p.value
                }
                #else {
                    #    male_p.value <- 10
                    #}
                if (!is.null(phenTestResult$model.output$female)){
                    RROutput <- phenTestResult$model.output$female
                    if (RROutput[1] < phenotypeThreshold && RROutput[3] >= phenotypeThreshold){
                        direction_females <- "Low"
                        female_p.value <- RROutput[1]
                    }
                    else if (RROutput[1] >= phenotypeThreshold && RROutput[3] < phenotypeThreshold){
                        direction_females <- "High"
                        female_p.value <- RROutput[3]
                    }
                    #female_p.value <- result$model.output$female$p.value
                }
                #else {
                    #    female_p.value <- 10
                    #}
                #all_p.value <- result$model.output$all$p.value
                
                ChangeClassification <- paste("Not significant")
                if (phenTestResult$numberSexes==2) {
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