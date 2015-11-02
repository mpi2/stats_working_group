## Copyright © 2012-2015 EMBL - European Bioinformatics Institute
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
## LogisticFramework.R contains LRDataset, startLRModel, finalLRModel, modelFormulaLR
## and parserOutputSummaryLR functions
##------------------------------------------------------------------------------
## LRDataset prepares dataset for the LR framework - maps values to 0/1
## Abnormal - 1 (will be modeled), Normal - 0 
LRDataset <- function(phenList=NULL, depVariable=NULL, abnormalValues=c("abnormal","Abnormal","TRUE","deviant"),
        outputMessages=TRUE)
{
    stop_message <- ""
    ## START CHECK ARGUMENTS   
    # 1
    if (is.null(phenList) || !is(phenList,"PhenList")){
        stop_message <- paste(stop_message,"Error:\nPlease create and specify PhenList object.\n",sep="")
    }
    # 2
    if (is.null(depVariable)){
        stop_message <- paste(stop_message,"Error:\nPlease specify dependent variable, ",
                "for example: depVariable='Lean.Mass'.\n",sep="")
    }
    if (nchar(stop_message)==0) {
        # create small data frame with values to analyse
        columnOfInterest <- getColumn(phenList,depVariable)
        # Presence
        if (is.null(columnOfInterest))
            stop_message <- paste("Error:\nDependent variable column '",
                depVariable,"' is missed in the dataset.\n",sep="")
    }    
    ## STOP CHECK ARGUMENTS   
    # If problems are deteckted
    if (nchar(stop_message)>0){
        if (outputMessages)  { 
            message(stop_message)
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt))      
            stop()
        }
        else {
            stop(stop_message)
        }
    }
    
    # Creates a subset with analysable variable
    if (nchar(stop_message)==0) {
        oldvalues <- levels(factor(na.omit(columnOfInterest)))
        oldvaluesAbnormal <- oldvalues[which(oldvalues %in% abnormalValues)]
        oldvaluesNormal <- oldvalues[which(!(oldvalues %in% abnormalValues))]
        newvalues <- c(rep(1,length(oldvaluesAbnormal)),rep(0,length(oldvaluesNormal)))
        map <- setNames(newvalues,c(oldvaluesAbnormal,oldvaluesNormal))
        mappedvalues <- sapply(columnOfInterest, function(x) map[as.character(x)])
        names(mappedvalues)<-NULL
        columnOfInterest <- as.integer(mappedvalues)
        phenList@datasetPL[,c(depVariable)] <- columnOfInterest
        return(phenList)
    }  
}
##------------------------------------------------------------------------------
## LRstartModel creates start model and modifies it after testing of different
## hypothesis (the model effects).
## The model effects are:
## - batch effect (random effect significance)
## Fixed effects:
## - interaction effect (genotype by sex interaction significance)
## - sex effect (sex significance)
##
## As the result PhenTestResult object that contains calculated or user
## defined model effects and MM start model is created.
##http://idiom.ucsd.edu/~rlevy/lign251/fall2007/lecture_15.pdf 
## States for nested mixed effect models liklihood ratio can be used to compare model
##provided only differ in fixed effects structure

#debug(startLRModel_SD)

startLRModel_SD <- function(phenList, depVariable, outputMessages=TRUE, pThreshold=0.05)

{
    x <- getDataset(phenList)
    numberofsexes <- length(levels(x$Sex))
    keep_batch <- NA
    keep_interaction <- NA
    keep_sex <- NA
    keep_equalvar <- NA  
    keep_weight <- NA  # weight not included as a fixed effect      
    equation <- "withoutWeight"   # weight not included as a fixed effect    
    
	
    #START OF tryCatch    
    finalResult <- tryCatch({                
                ## We're using the same framework as for LR and MM.  typically here would be the model optmiisation eg batch/heterogenous variance.  In this pipeline we are not completing any optimisation but we are using the framework.     
                                                      
                ## Need to obtain final SD model 
				#x <- getDataset(phenList)  
				KOdataOnly <- subset(x,x$Genotype==refGenotype(phenList))  #SD assessment only KO data
				L_model_SexOnlyModel = modelFormulaLR_SD(depVariable, sexIncluded=TRUE, dimorphismIncluded=FALSE, genotypeIncluded=FALSE)
                finalSDModel <- do.call("logistf", args = list(L_model_SexOnlyModel, data =KOdataOnly , na.action="na.omit"))
                
				#LR modelling output 
                linearRegressionOutput <- list(model.output=finalSDModel,
                        equation = equation,
                        model.effect.batch=keep_batch,
                        model.effect.interaction=FALSE,
                        model.output.interaction=NA,
                        model.effect.sex=keep_sex,
                        model.effect.weight=keep_weight,
                        model.effect.variance=NA,
                        numberSexes=numberofsexes,
                        pThreshold=pThreshold,
                        model.formula.genotype=L_model_SexOnlyModel)
                linearRegressionOutput$FET <- FisherExactTest(phenList,depVariable) 
                
                 finalResult <- new("PhenTestResult",
                    analysedDataset=x,
                    depVariable=depVariable,
                    refGenotype=refGenotype(phenList),
                    testGenotype=testGenotype(phenList),
                    method="SD_categorical",
                    transformationRequired = FALSE,
                    lambdaValue = integer(0),
                    scaleShift = integer(0),
                    analysisResults=linearRegressionOutput) 
            },
            #END OF tryCatch    
            error=function(error_mes) {
                finalResult <- NULL                
                stop_message <- paste("Error:\nCan't fit the model ",
                        format(formula_withOutBatch),". Try FE method.\n",sep="")                  
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
    )                    
    return(finalResult)    
}
##------------------------------------------------------------------------------
## Function to build start formula for Logistic framework 
modelFormulaLR_SD <- function(depVariable, genotypeIncluded, sexIncluded, dimorphismIncluded)
{         
                
                    if (!sexIncluded  && !dimorphismIncluded && !genotypeIncluded){
						LR_model_formula =  as.formula(paste(depVariable, "~", "1", sep=""))
                    } else if(sexIncluded && dimorphismIncluded && genotypeIncluded){
						LR_model_formula =   as.formula(paste(depVariable, "~",paste("Genotype", "Sex",  "Genotype*Sex", sep= "+")))
                    } else if(sexIncluded  && !dimorphismIncluded && !genotypeIncluded){
						LR_model_formula =   as.formula(paste(depVariable, "~","Sex", sep= ""))
					} else if(sexIncluded && !dimorphismIncluded && genotypeIncluded){
						LR_model_formula =   as.formula(paste(depVariable, "~",paste("Genotype", "Sex",  sep= "+")))
					}
                
    return(LR_model_formula)
}

##------------------------------------------------------------------------------
## Works with PhenTestResult object created by testDataset function.
## Builds null and uses final model to assess genotype effect
## final model results are then captured and stored.  
finalLRModel_SD <- function(phenTestResult, outputMessages=TRUE)
{
    stop_message <- ""    
    ## START Check PhenTestResult object
    if(is(phenTestResult,"PhenTestResult")) {
		
		#preparation of various datasets
		x <- analysedDataset(phenTestResult)
		WTdataOnly <- subset(x,x$Genotype==testGenotype(phenTestResult))
		KOdataOnly <- subset(x,x$Genotype==refGenotype(phenTestResult))
		
		result <- phenTestResult  
        linearRegressionOutput <- analysisResults(result)
        
		depVariable <- getVariable(result)
        numberofsexes <- linearRegressionOutput$numberSexes
        keep_sex <- linearRegressionOutput$model.effect.sex
        keep_interaction <- FALSE
        keep_batch <- linearRegressionOutput$model.effect.batch        
       
		## Stop function if there are no datasets to work with
        if(is.null(x))
            stop_message <- "Error:\nPlease create a PhenList object first and run function 'testDataset'.\n"      
        ## Stop function if there are no enough input parameters
        if ( is.null(depVariable) || is.null(keep_batch)
                || is.null(keep_sex) || is.null(keep_interaction))
            stop_message <- "Error:\nPlease run function 'testDataset' first.\n"
    }
    else{
        stop_message <- "Error:\nPlease create a PhenTestResult object first.\n"
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
    ## END Checks and stop messages    
    result <- tryCatch( {   
                ## Test: genotype groups association with dependent variable
                          
				#various formula's needed to build the various models
				L_model_FullModel = modelFormulaLR_SD(depVariable, sexIncluded=TRUE, dimorphismIncluded=TRUE, genotypeIncluded=TRUE)
				L_model_SexOnlyModel = modelFormulaLR_SD(depVariable, sexIncluded=TRUE, dimorphismIncluded=FALSE, genotypeIncluded=FALSE)
				L_model_FullwithoutInteraction = modelFormulaLR_SD(depVariable, sexIncluded=TRUE, dimorphismIncluded=FALSE, genotypeIncluded=TRUE)
				      
			 	#Build various models 
				#all data:
                Full_Model <- do.call("logistf", args = list(L_model_FullModel, data = x, na.action="na.omit")) 
				SexOnly_Model<- do.call("logistf", args = list(L_model_SexOnlyModel, data = x, na.action="na.omit")) 
				#WT data:
				SexOnly_WTonly_Model<- do.call("logistf", args = list(L_model_SexOnlyModel, data = WTdataOnly, na.action="na.omit")) 
				#KO data:
				SexOnly_KOonly_Model<- do.call("logistf", args = list(L_model_SexOnlyModel, data = KOdataOnly, na.action="na.omit")) 
				
				#TestingImpact
				#Sex
				ImpactSex_p.value=logistftest(SexOnly_WTonly_Model)$prob     #can't build an intercept only model with logistif so this the route to compare against intercept only model
				#Genotype
				genotypeTest_p.value = anova(Full_Model, SexOnly_Model)$pval  
				#SD
				SD_p.value=logistftest(SexOnly_KOonly_Model)$prob   #can't build an intercept only model with logistif so this the route to compare against intercept only model
								
				#storing the results			
               	#those associated with testing of SD
				linearRegressionOutput$SD.model.output <-  SexOnly_KOonly_Model 
				linearRegressionOutput$SD.model.formula <- L_model_SexOnlyModel
				linearRegressionOutput$SD.model.output.pVal <- SD_p.value
				#those associated with testing of Sex
				linearRegressionOutput$Sex.model.output <-  SexOnly_WTonly_Model
				linearRegressionOutput$Sex.model.formula <- L_model_SexOnlyModel
				linearRegressionOutput$Sex.model.output.pVal <- ImpactSex_p.value
				#those associated with testing of Sex
				linearRegressionOutput$Genotype.model.output <-  Full_Model
				linearRegressionOutput$Genotype.model.formula <- L_model_FullModel
				linearRegressionOutput$Genotype.model.output.pVal <- genotypeTest_p.value

								
				# Create modeloutput and choose output depending on model
                linearRegressionOutput$model.output.summary  <- parserOutputSummaryLR_SD(linearRegressionOutput)  
                
				
				#store all these results approproately in the results object 
				result@analysisResults <- linearRegressionOutput 
               
				#Assign LR quality of fit
                result@analysisResults$model.output.quality <- testFinalLRModel(result)
                result
            },            
            # End of tryCatch statement - if fails try to suggest something useful for the user
            error=function(error_mes) {
                result <- NULL
                stop_message <- paste("Error:\nCan't fit the model ", format(result$L_model_FullModel),
                        ". Doh!.\n",sep="")
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
    )        
    return(result)
}
##------------------------------------------------------------------------------
## Creates model output summary and returns it in readable vector format
parserOutputSummaryLR_SD<-function(linearRegressionOutput)
{
	result <- linearRegressionOutput
	modeloutput_summary <- summary(result$model.SDmodel_output)
	
	# Estimates of coefficients for the interaction model
    coefficients <- format(linearRegressionOutput$SD.model.output$coefficients,scientific=FALSE)
    # Standard errors for the interaction model
    error_estimates=format(sqrt(diag(vcov(linearRegressionOutput$SD.model.output ))),scientific=FALSE)
    # p-values for the interaction model
    probs <- format(linearRegressionOutput$SD.model.output$prob,scientific=FALSE)        
    
	# Set all values to NA initially prior to selecting those relevant to the model   # matching output product to that used in the continuous pipeline hence the empty terms
    genotype_estimate <- NA
    genotype_estimate_SE <- NA
    genotype_p_value <- NA    
    sex_estimate <- NA
    sex_estimate_SE <- NA
    sex_p_value <- NA    
    intercept_estimate <- NA
    intercept_estimate_SE <- NA
    weight_estimate <- NA
    weight_estimate_SE <- NA
    weight_p_value <- NA    
    sex_FvKO_estimate <- NA
    sex_FvKO_SE <- NA
    sex_FvKO_p_value <- NA
    sex_MvKO_estimate <- NA
    sex_MvKO_SE <- NA
    sex_MvKO_p_value <- NA    
    
	#Capture model values
	intercept_estimate <- coefficients[1]
    intercept_estimate_SE <- error_estimates[1]    
    sex_estimate <- coefficients[2]
    sex_estimate_SE <- error_estimates[2]
	sex_p_value <- probs[2]    
        
    output <- c(genotype_estimate, genotype_estimate_SE, genotype_p_value, sex_estimate, sex_estimate_SE,  sex_p_value, 
            "NA", "NA", "NA", 
            intercept_estimate, intercept_estimate_SE, 
            sex_FvKO_estimate, sex_FvKO_SE, sex_FvKO_p_value,  
            sex_MvKO_estimate, sex_MvKO_SE, sex_MvKO_p_value)    
    names(output) <- c("genotype_estimate", "genotype_estimate_SE", 
            "genotype_p_value", 
            "sex_estimate", "sex_estimate_SE", "sex_p_value", 
            "weight_estimate", "weight_estimate_SE", "weight_p_value", 
            "intercept_estimate", "intercept_estimate_SE", 
            "sex_FvKO_estimate", "sex_FvKO_SE", "sex_FvKO_p_value", 
            "sex_MvKO_estimate", "sex_MvKO_SE", "sex_MvKO_p_value")
    return(output) 
}
##------------------------------------------------------------------------------