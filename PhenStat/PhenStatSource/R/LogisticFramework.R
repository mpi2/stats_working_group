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
## LogisticFramework.R contains startModel, finalModel, modelFormula
## and parserOutputSummary functions
##------------------------------------------------------------------------------
## startModel creates start model and modifies it after testing of different
## hypothesis (the model effects).
## The model effects are:
## - batch effect (random effect significance)
## Fixed effects:
## - interaction effect (genotype by sex interaction significance)
## - sex effect (sex significance)
##
## As the result PhenTestResult object that contains calculated or user
## defined model effects and MM start model is created.


###Using Firth penalized logistic regression designed for low n with spearation issues.  Note cannot deal with random effect models so am testing for batch but then proceeding with standard model. 


LR_Model <- function(phenList, depVariable, 	outputMessages=TRUE, pThreshold=0.0, baselineLevel)
{
	x <- phenList$dataset
	numberofsexes <- length(levels(x$Sex))
			
	## Start model formula: genotype, sex and genotype*sex interaction included.  Note if only one sex then sex and interaction omitted 
	##very different to MM as random effect built in here into formula. 
	formula_withBatch=modelFormula_LF(numberofsexes, depVariable, sexIncluded=TRUE, dimorphismIncluded=TRUE, IncludeBatch="Yes")
	formula_withOutBatch=modelFormula_LF(numberofsexes, depVariable, sexIncluded=TRUE, dimorphismIncluded=TRUE, IncludeBatch="No")		
	
	#setlevels
	 x$Genotype=relevel(x$Genotype, ref="+/+")
	 x[ , depVariable]=relevel(x[ ,depVariable], ref=baselineLevel)
	
	#START OF tryCatch    
	finalResult <- tryCatch({
				
				
				##Goal of this section is to assess whether batch is significant or not in explaining variation
				if ('Batch' %in% colnames(x)){
					## LR fit of model formula with no random effects
					## Model 1A (model_withoutbatch)
					L_model_withoutbatch <- do.call("glm",args=list(formula_withOutBatch, x, na.action="na.omit", family=binomial()  ))
						
					## MM fit of model formula (with random effects)
					## Model 1 (model_withBatch)
					L_model_withBatch <- do.call("glmer", args = list(formula_withBatch, data = x, na.action="na.omit",family=binomial() ))	
										
					## Test: the random effects associated with batch intercepts can be ommited from model
					## Hypothesis 1
					## Null Hypothesis: variance of batch = 0
					## Alternative Hypothesis: variance of batch > 0	
					##If p value below threshold then you reject null and accept alternative that batch is significant in explaining variation in the model
					##Based on method shown here http://www.simonqueenborough.com/R/specialist/mixed-models.html
					
					#p.value.batch = anova(L_model_withBatch, L_model_withoutbatch, test = "Chisq")  NOT suitable for comparing mixed against standard model
					p.value.batch <- pchisq(-2*(logLik(L_model_withBatch)-logLik(L_model_withoutbatch)), 1, lower=FALSE)[1] 
					keep_batch <- p.value.batch<pThreshold
					##within this framework we are not testing the concept that the variance might depend on the genotype group
					keep_equalvar <- NA	
				}
				else {
					## No Batch effects
					keep_batch <- FALSE
					keep_equalvar <- NA	
				}
				
				## Goal of this section is to tests for significance of fixed effects by comparing models with and without fixed effects included		
				if(numberofsexes==2){
					
									
						#test sexual dimorphism
						L_model_withoutbatch <- do.call("logistf",	args=list(formula_withOutBatch, x, na.action="na.omit", family=binomial()))
						formula_withoutBatch_noSD=modelFormula_LF(numberofsexes, depVariable, sexIncluded=TRUE, dimorphismIncluded=FALSE, IncludeBatch="No")
						L_model_withoutBatch_noSD <- do.call("logistf", args = list(formula_withoutBatch_noSD, data = x, na.action="na.omit" ))
						
						interactionTest = anova(L_model_withoutbatch, L_model_withoutBatch_noSD)$pval
						keep_interaction <- interactionTest<pThreshold
						
						#test sex inclusion in model
						formula_withoutBatch_noSex=modelFormula_LF(numberofsexes, depVariable, sexIncluded=FALSE, dimorphismIncluded=FALSE, IncludeBatch="No")
						L_model_withoutBatch_noSex <- do.call("logistf", args = list(formula_withoutBatch_noSex, data = x, na.action="na.omit"))
						
						sexTest = anova(L_model_withoutbatch, L_model_withoutBatch_noSex)$pval
						keep_sex <- sexTest<pThreshold
						keep_weight <- NA														
																	
				}
				else {
					keep_sex <- FALSE
					keep_interaction <- FALSE
					interactionTest <- NA
					keep_weight <- NA
				}
				
				if (outputMessages)
					message(paste("Information:\nCalculated values for model effects are: ",
									"keepBatch=",keep_batch,
									", keepSex=",keep_sex,
									", keepInteraction=",keep_interaction,".\n",sep=""))
				
				
				## Need to obtain final model 
				##step one build formula and then fit the model through glm or glmer route depending on whether batch is significant
				
									
				FinalFormula=genotype_modelFormula_LF(numberofsexes, depVariable, sexIncluded=keep_sex, dimorphismIncluded=keep_interaction, IncludeBatch="No")	
				finalModel <- do.call("logistf", args = list(FinalFormula, data = x, na.action="na.omit"))
									
						
						
				finalResult <- new("PhenTestResult",list(
								model.dataset=x,
								model.output=finalModel,
								depVariable=depVariable,
								equation=NA,
								method="LR",
								model.effect.batch=keep_batch,
								model.effect.variance=NA,
								model.effect.interaction=keep_interaction,
								model.output.interaction=interactionTest,
								model.effect.sex=keep_sex,
								model.effect.weight=NA,
								numberSexes=numberofsexes,
								pThreshold=pThreshold,
								model.formula.genotype=FinalFormula,   
								model.output.averageRefGenotype=NA))
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

############################################function written by Natasha to build start formula for Logistic framework which includes batch specified in the formula to allow the various models to be compared

modelFormula_LF <- function(numberofsexes, depVariable, sexIncluded, dimorphismIncluded, IncludeBatch){
	
	LR_model_formula <- switch(IncludeBatch,
			Yes = {
				## mixed logistic model framework
				if(numberofsexes==2){
					if (!sexIncluded  && !dimorphismIncluded){
						as.formula(paste(depVariable, "~", paste("Genotype", "(1|Batch)", sep="+")))
					} else if((sexIncluded && dimorphismIncluded)|(!sexIncluded && dimorphismIncluded)){
						as.formula(paste(depVariable, "~",paste("Genotype", "Sex",  "Genotype*Sex", "(1|Batch)", sep= "+")))
					} else if(sexIncluded  && !dimorphismIncluded){
						as.formula(paste(depVariable, "~",paste("Genotype","Sex","(1|Batch)", sep= "+")))
					}
					
				}else{
					as.formula(paste(depVariable, "~", paste("Genotype", "(1|Batch)", sep="+")))
					}
			},
			No = {
				## standard logistic model framework
				if(numberofsexes==2){
					if (!sexIncluded  && !dimorphismIncluded){
						as.formula(paste(depVariable, "~", "Genotype", sep=""))
					} else if((sexIncluded && dimorphismIncluded)|(!sexIncluded && dimorphismIncluded)){
						as.formula(paste(depVariable, "~",paste("Genotype", "Sex",  "Genotype*Sex", sep= "+")))
					} else if(sexIncluded  && !dimorphismIncluded){
						as.formula(paste(depVariable, "~",paste("Genotype","Sex", sep= "+")))
					}
				}else{
					as.formula(paste(depVariable, "~","Genotype", sep=" "))
				}
			}
	)
	return(LR_model_formula)
}


null_modelFormula_LF <- function(numberofsexes, depVariable, sexIncluded, dimorphismIncluded, IncludeBatch){
	
	null_LR_model_formula <- switch(IncludeBatch,
			Yes = {
				## mixed logistic model framework
				if(numberofsexes==2){
					if (!sexIncluded  && !dimorphismIncluded){
						as.formula(paste(depVariable, "~", paste(1, "(1|Batch)", sep="+")))
					} else if((sexIncluded && dimorphismIncluded)|(!sexIncluded && dimorphismIncluded)|(sexIncluded  && !dimorphismIncluded)){
						as.formula(paste(depVariable, "~",paste("Sex",  "(1|Batch)", sep= "+")))
					}
					
				}else{
					as.formula(paste(depVariable, "~", paste(1, "(1|Batch)", sep="+")))
				}
			},
			No = {
				## standard logistic model framework
				if(numberofsexes==2){
					if (!sexIncluded  && !dimorphismIncluded){
						as.formula(paste(depVariable, "~", 1, sep=""))
					} else if((sexIncluded && dimorphismIncluded)|(!sexIncluded && dimorphismIncluded)|(sexIncluded  && !dimorphismIncluded)){
						as.formula(paste(depVariable, "~",paste( "Sex", sep= "+")))
					}
				}else{
					as.formula(paste(depVariable, "~", "1", sep=" "))
				}
			}
			
	)
	return(null_LR_model_formula)
}

#Second cycle needed as specification of model changes for the final fitting in the presence of sexual dimorphism
genotype_modelFormula_LF <- function(numberofsexes, depVariable, sexIncluded, dimorphismIncluded, IncludeBatch){
	
	LR_model_formula <- switch(IncludeBatch,
			Yes = {
				## mixed logistic model framework
				if(numberofsexes==2){
					if (!sexIncluded  && !dimorphismIncluded){
						as.formula(paste(depVariable, "~", paste("Genotype", "(1|Batch)", sep="+")))
					} else if((sexIncluded && dimorphismIncluded)|(!sexIncluded && dimorphismIncluded)){
						as.formula(paste(depVariable, "~",paste("Sex",  "Genotype:Sex", "(1|Batch)", sep= "+")))
					} else if(sexIncluded  && !dimorphismIncluded){
						as.formula(paste(depVariable, "~",paste("Genotype","Sex","(1|Batch)", sep= "+")))
					}
					
				}else{
					as.formula(paste(depVariable, "~", paste("Genotype", "(1|Batch)", sep="+")))
				}
			},
			No = {
				## standard logistic model framework
				if(numberofsexes==2){
					if (!sexIncluded  && !dimorphismIncluded){
						as.formula(paste(depVariable, "~", "Genotype", sep=""))
					} else if((sexIncluded && dimorphismIncluded)|(!sexIncluded && dimorphismIncluded)){
						as.formula(paste(depVariable, "~",paste("Sex",  "Genotype:Sex", sep= "+")))
					} else if(sexIncluded  && !dimorphismIncluded){
						as.formula(paste(depVariable, "~",paste("Genotype","Sex", sep= "+")))
					}
				}else{
					as.formula(paste(depVariable, "~","Genotype", sep=" "))
				}
			}
	)
	return(LR_model_formula)
}

###------------------------------------------------------------------------------
## Works with PhenTestResult object created by testDataset function.
## Builds null and uses final model to assess genotype effect
## final model results are then captured and stored.  

queryFinalModel <- function(phenTestResult, outputMessages=TRUE, baselineLevel)
{
	## Checks and stop messages
	stop_message <- ""
	
	## Check PhenTestResult object
	if(is(phenTestResult,"PhenTestResult")) {
		result <- phenTestResult
		x <- result$model.dataset
		depVariable <- result$depVariable
		x$Genotype=relevel(x$Genotype, ref="+/+")
		x[ , depVariable]=relevel(x[ ,depVariable], ref=baselineLevel)
		
		
		
		equation <- NA
		keep_weight <- NA
		keep_sex <- result$model.effect.sex
		keep_interaction <- result$model.effect.interaction
		keep_batch <- result$model.effect.batch
		keep_equalvar <- NA
		
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
	
	numberofsexes <- result$numberSexes
	
	finalResult <- tryCatch( {   
				## Test: genotype groups association with dependent variable
				## Null Hypothesis: genotypes are not associated with dependent variable
				## Alternative Hypothesis: genotypes are associated with dependent  variable
				#setlevels
							
				
				#build null model
					model_null.formula = null_modelFormula_LF (numberofsexes, depVariable, sexIncluded=keep_sex, dimorphismIncluded=keep_interaction, IncludeBatch="No")					
					nullModel <- do.call("logistf", args = list(model_null.formula, data = x, na.action="na.omit"))
					#build genotype model
					model_genotype.formula = genotype_modelFormula_LF (numberofsexes, depVariable, sexIncluded=keep_sex, dimorphismIncluded=keep_interaction, IncludeBatch="No")
					genotypeModel <- do.call("logistf", args = list(model_genotype.formula, data = x, na.action="na.omit" ))						
					#compare with genotype model
					#genotypeTest_p.value <- pchisq((deviance(genotypeModel)-deviance(nullModel)), 1, lower=FALSE)  #problem with the estimate of df
					genotypeTest_p.value = anova(genotypeModel, nullModel)$pval 
					
				
								
				## Store the results
				result$model.genotype <-genotypeModel
				result$model.formula.genotype <- model_genotype.formula
				result$model.null <-nullModel
				result$model.formula.null <- model_null.formula
				result$model.output.genotype.nulltest.pVal <- genotypeTest_p.value 
				result$model.effect.variance <- NA
				
				## Assign MM quality of fit
				result$model.output.quality <- testFinal_LR_Model(phenTestResult)   
				
				## Parse modeloutput and choose output depending on model
				result$model.output.summary <- parserOutputSummary(result)
				
			},
			
			# End of tryCatch statement - if fails try to suggest something useful for the user
			error=function(error_mes) {
				
				finalResult <- NULL
				stop_message <- paste("Error:\nCan't fit the model ",	format(result$model_genotype.formula),". Try FE Framework.\n",sep="")
				
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
	finalResult <- result
	
	return(finalResult)
}







##------------------------------------------------------------------------------
## Parses model output summary and returns in readable vector format
parserOutputSummary<-function(phenTestResult)
{
	result <- phenTestResult
	modeloutput_summary <- summary(result$model.output)
	
	#set all values to NA initially prior to selecting those relevant to the model 
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
	
	#pull out the information that is used to drive decision tree in where the information resides
	keep_sex <- result$model.effect.sex
	keep_interaction <- result$model.effect.interaction
	keep_batch <- result$model.effect.batch
	no_of_sexes <- result$numberSexes
	
	#decision tree to pull the information depending on the final model fitted
	#note position is not dependent on whether a mixed or standard logisitic model is fitted.
	if(no_of_sexes==1){
		
		intercept_estimate = modeloutput_summary$coefficients[1]
		intercept_estimate_SE 	= modeloutput_summary$coefficients[3]
		
		genotype_estimate = modeloutput_summary$coefficients[2]
		genotype_estimate_SE = modeloutput_summary$coefficients[4]
		genotype_p_value =  modeloutput_summary$coefficients[8]
		
	}else if(keep_interaction==TRUE){
					
	
		intercept_estimate = modeloutput_summary$coefficients[1]
		intercept_estimate_SE 	= modeloutput_summary$coefficients[5]
					
		sex_estimate=modeloutput_summary$coefficients[2]
		sex_estimate_SE=modeloutput_summary$coefficients[6]
		sex_p_value= modeloutput_summary$coefficients[14]
		
		sex_FvKO_estimate= modeloutput_summary$coefficients[3]
		sex_FvKO_SE=modeloutput_summary$coefficients[7]
		sex_FvKO_p_value=modeloutput_summary$coefficients[15]
		sex_MvKO_estimate=modeloutput_summary$coefficients[4]
		sex_MvKO_SE=modeloutput_summary$coefficients[8]
		sex_MvKO_p_value=modeloutput_summary$coefficients[16]		
	
	}else{
		
		intercept_estimate = modeloutput_summary$coefficients[1]
		intercept_estimate_SE 	= modeloutput_summary$coefficients[4]
		
		sex_estimate=modeloutput_summary$coefficients[3]
		sex_estimate_SE=modeloutput_summary$coefficients[6]
		sex_p_value= modeloutput_summary$coefficients[12]
		
		genotype_estimate = modeloutput_summary$coefficients[2]
		genotype_estimate_SE = modeloutput_summary$coefficients[5]
		genotype_p_value =  modeloutput_summary$coefficients[11]
		
	}	
		
	output <- c(genotype_estimate, genotype_estimate_SE,  genotype_p_value,
			sex_estimate, sex_estimate_SE,  sex_p_value, 
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
