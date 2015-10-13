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
        if (method(phenTestResult) %in% c("MM","TF","LR", "SD_continuous")){
            equation <- analysisResults$equation
            keep_weight <- analysisResults$model.effect.weight
            keep_sex <- analysisResults$model.effect.sex
            keep_interaction <- analysisResults$model.effect.interaction
            keep_batch <- analysisResults$model.effect.batch
            keep_equalvar <- analysisResults$model.effect.variance
            model.output <- analysisResults$model.SDmodel_output
            
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
        # SD
        if (method(phenTestResult) %in% c("SD_continuous" )) {
            if (userMode=="summaryOutput") {
                if(is.na(analysisResults$model.output.SDtest.pVal)==TRUE){
                    ChangeClassification <- NA
                }else if(analysisResults$model.output.SDtest.pVal>phenotypeThreshold){
                    
                        ChangeClassification <- paste("With phenotype threshold value",
                                phenotypeThreshold,"- no significant change")
                    
				}else  if(analysisResults$model.output.summary["sex_FvKO_p_value"]>=0.05 
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
            
            #vectorOutput
            else {
                
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
        
        
        return(ChangeClassification)
        
   	 }
	}

}

##------------------------------------------------------------------------------