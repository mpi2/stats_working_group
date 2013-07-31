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
# testDataset.R contains testDataset and modelFormula functions and constructs the PhenTestResult 
# class object
#-----------------------------------------------------------------------------------
testDataset <- function(phenList, depVariable, equation="withWeight", outputMessages=TRUE, pThreshold=0.05, method="MM", keepList=NULL)
# Performs the following checks for dependent variable:
# 1. "depVariable" column should present in the dataset
# 2. "depVariable" should be numeric for Mixed Model (MM) framework, otherwise recommends Fisher Exact Test (FE)
# 3. Each one genotype level should have more than one "depVariable" level
# (variability) for MM framework, otherwise recommends FE framework. 

# Performs other consistency checks. For instance, when "depVariable" column contains cathegorical values, 
# but defined method is MM redefines method to FE.

# MM framework: 
# For MM framework creates start model and modifies it after testing of different hypothesis (the model effects).
# The model effects are: 
# - batch effect (random effect significance), 
# - variance effect (TRUE if residual variances for genotype groups are homogeneous and FALSE if they are heterogeneous), 
# - interaction effect (genotype by gender interaction significance), 
# - gender effect (gender significance), 
# - weigth effect (weigth significance).

# If user would like to assign other TRUE/FALSE values to the effects of the model 
# then he or she has to define keepList argument which is vector of TRUE/FALSE values.

# If user has defined model effects (keepList argument) then function prints out calculated and user defined effects 
# (only when outputMessages argument is set to TRUE), checks user defined effects for consistency 
# (for instance, if there are no "Weight" column in the dataset then weigth effect can't be assigned to TRUE, etc.)
# and modifies start model according to user defined effects. 

# As the result PhenTestResult object that contains calculated or user defined model effects and MM start model is created. 

# FE framework:  
# For Fisher Exact Test framework there are no special checks except "depVariable" checks mentioned above. 
# "buildFisherExactTest" function is called for count matrix (matrices) creation and Fihser tests performance. 

# As the result PhenTestResult object that contains Fisher tests output is created.   

{
    require(nlme)
    
    stop_message <- ""
    
    # Checks and stop messages
    
    if(is(phenList,"PhenList")) {
        x <- phenList$dataset    
        
    } else {
        stop_message <- "Error:\nPlease create a PhenList object first.\n"
    }    
    
    if (!(depVariable %in% colnames(x)))
    stop_message <- paste("Error:\nDependent variable column '",depVariable,"' is missed in the dataset.\n",sep="")
    
    
    if (!(equation %in% c("withWeight","withoutWeight")) && method=="MM")
    stop_message <- "Error:\nPlease define equation you would like to use from the following options: 'withWeight', 'withoutWeight'\n."
    
    if (!('Weight' %in% colnames(x)) && equation=="withWeight" && method=="MM"){
        if (outputMessages)
        message("Warning:\nWeight column is missed in the dataset. Equation 'withWeight' can't be used and has been replaced to 'withoutWeight'.")
        equation="withoutWeight"
    }
    
    # Test: depVariable is continuous variable
    
    columnOfInterest <- x[,c(depVariable)]
    
    if(is.numeric(columnOfInterest)){
        if ((length(unique(columnOfInterest))/length(columnOfInterest)<0.05) && outputMessages && method=="MM") 
        message(paste("Warning: Dependent variable '",depVariable,"' is numeric but seemed to be categorical because there is little variation. Fisher Exact Test can be better way to do the analysis than Mixed Models.\n",sep="")) 
    }
    else if (method=="MM"){
        method="FE"
        if (outputMessages)
        message(paste("Warning:\nDependent variable '",depVariable,"' is not numeric. Fisher Exact Test will be used for the analysis of this dependent variable.\n",sep=""))
    }
    # Test: depVariable variablity in Genotypes (require at least 2 levels)
    
    Genotype_levels=levels(x$Genotype)
    Gender_levels=levels(x$Gender)  
    
    for (i in 1:length(Genotype_levels)){
        GenotypeSubset <- subset(x, x$Genotype==Genotype_levels[i])
        for (j in 1:length(Gender_levels)){           
            GenotypeGenderSubset <- subset(GenotypeSubset, GenotypeSubset$Gender==Gender_levels[j]) 
            columnOfInterest <- GenotypeGenderSubset[,c(depVariable)]
            if (length(unique(columnOfInterest))==1)
            stop_message <- paste("Error:\nInsufficient variability in the dependent variable '",depVariable,"' for genotype/gender combinations to allow the application of Mixed Model or Fisher Exact test framework.\n",sep="")
        }                
    }  
    
    # Dealing with provided significance values
    if (!is.null(keepList)){
        # Stop function if there are no enough needed input parameters
        if (length(keepList)!=5) 
        stop_message <- "Error:\nPlease define the values for 'keepList' list, where for each effect/part of the model TRUE/FALSE value defines to keep it in the model or not: 
        'keepList=c(keepBatch,keepVariance,keepWeight,keepGender,keepInteraction)'.\n"
        
    }
    
    
    if (nchar(stop_message)>1){
        if (outputMessages)   
        message(stop_message)
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))      
        stop()
    }
    
    # END Checks and stop messages
    
    if (outputMessages)
        message(paste("Information:\nDependent variable: '",depVariable,"'.\n",sep="")) 
    
    # Mixed Models framework
    if (method=="MM")    { 
        
        if (outputMessages)
            message(paste("Information:\nMethod: Mixed Model framework.\n",sep="")) 
        
        if (!is.null(keepList)){
            
            if (!('Weight' %in% colnames(x)) && keep_weight){
                if (outputMessages)
                    message("Warning:\nWeight column is missed in the dataset. 'keepWeight' is set to FALSE.")
                user_keep_weight=FALSE
            }
            
            if (!('Batch' %in% colnames(x)) && (keep_batch || keep_equalvar)){
                if (outputMessages)
                    message("Warning:\nBatch column is missed in the dataset. 'keepBatch' and 'keepVariance' are set to FALSE.")
                user_keep_batch=FALSE
                user_keep_equalvar=FALSE
            }
        
            # User's values for effects
            user_keep_weight <- keepList[3]
            user_keep_gender <- keepList[4]
            user_keep_interaction <- keepList[5]
            user_keep_batch <- keepList[1]
            user_keep_equalvar <- keepList[2]
            
            if (outputMessages)
                message(paste("Information:\nUser's values for model effects are: keepBatch=",user_keep_batch,
                            ", keepVariance=",user_keep_equalvar,
                            ", keepWeight=",user_keep_weight,
                            ", keepGender=",user_keep_gender, 
                            ", keepInteraction=",user_keep_interaction,".\n",sep=""))     
        }
        
        numberofgenders=length(levels(x$Gender))
        
        
        # Start model formula: homogenous residual variance, genotype and sex interaction included  
        model.formula = modelFormula(equation,numberofgenders, depVariable)
        
        if ('Batch' %in% colnames(x)){
            
            # MM fit of model formula (with random effects)
            # Model 1 (model_MM)
            model_MM = do.call("lme", args = list(model.formula, random=~1|Batch, data = x, na.action="na.omit", method="REML"))
            # GLS fit of model formula  (no random effects)
            # Model 1A (model_withoutbatch)
            model_withoutbatch <- do.call("gls", args=list(model.formula, x, na.action="na.omit"))
            # Test: the random effects associated with batch intercepts can be ommited from model
            # Hypothesis 1
            # Null Hypothesis: variance of batch = 0
            # Alternative Hypothesis: variance of batch > 0 
            # For the division by 2 explanations see p.80 of "Linear Mixed Models"...
            p.value.batch <-(anova(model_MM, model_withoutbatch)$p[2])/2
            # The result of the test for Hypothesis 1 will help to select the structure for random effects
            keep_batch= p.value.batch<pThreshold
            
            # MM fit of model formula with heterogeneous residual variances for genotype groups
            # Model 1 assumes homogeneous residual variances
            # Model 2 with heterogeneous residual variances 
            model_hetvariance= do.call("lme", args=list(model.formula, random=~1|Batch, x, weights=varIdent(form=~1|Genotype), na.action="na.omit", method="REML"))
            # Test: the variance of the residuals is the same (homogeneous) for all genotype groups
            # Hypothesis 2
            # Null Hypothesis: all residual variances are equal
            # Alternative Hypothesis: the residue variance is not equal    
            p.value.variance=(anova(model_MM, model_hetvariance)$p[2])
            # The result of the test for Hypothesis 2 will help to select a covariance structure for the residuals
            keep_equalvar= p.value.variance>pThreshold
        }
        else {
            # No Batch effects
            keep_batch = FALSE
            keep_equalvar = FALSE
        }
        
        
        # Model fit is selected according to test results    
        if(keep_batch && keep_equalvar){
            # Model 1
            model= model_MM
        }else if(keep_batch && !keep_equalvar){
            # Model 2
            model= model_hetvariance
        }else if(!keep_batch && keep_equalvar){
            # Model 1A
            model= model_withoutbatch
        }else if(!keep_batch && !keep_equalvar){
            # Modify model 1A to heterogeneous residual variances
            model= do.call("gls", args=list(model.formula, weights=varIdent(form=~1|Genotype), x, na.action="na.omit"))
        }
        
        
        # Tests for significance of fixed effects using TypeI F-test from anova functionality by using selected model
        anova_results = anova(model, type="marginal")$"p-value" < pThreshold
        if(numberofgenders==2){
            # Result of the test for gender significance (fixed effect 1.)
            keep_gender = anova_results[3]
            # Eq.2
            if (equation=="withWeight"){   
                # Result of the test for weight significance  (fixed effect 3.)        
                keep_weight = anova_results[4]
                # Result of the test for genotype by gender interaction significance (fixed effect 2.)
                keep_interaction = anova_results[5]
                
                # Technical results needed for the output
                # Interaction test results are kept for the output 
                interactionTest=anova(model, type="marginal")$"p-value"[5]           
            }
            # Eq.1
            else{
                # Result of the test for weight significance  (fixed effect 3.) 
                # It's FALSE since here the equation 1 is used - without weight effect
                keep_weight = FALSE
                # Result of the test for genotype by gender interaction significance (fixed effect 2.)
                keep_interaction = anova_results[4]
                # Interaction test results are kept for the output
                interactionTest=anova(model, type="marginal")$"p-value"[4]      
            } 
        }
        else {
            keep_gender = FALSE
            keep_interaction = FALSE 
            interactionTest = NA
            if (equation=="withWeight") keep_weight = anova_results[3]
            else keep_weight = FALSE
        }
        
        if (!keep_weight && equation=="withWeight") {
            equation="withoutWeight"
            if (outputMessages)
            message("Since weight effect is not significant the equation Eq.1 'withoutWeight' should be used instead.")
        }
        
        if (outputMessages)
            message(paste("Information:\nCalculated values for model effects are: keepBatch=",keep_batch,
                        ", keepVariance=",keep_equalvar,
                        ", keepWeight=",keep_weight,
                        ", keepGender=",keep_gender, 
                        ", keepInteraction=",keep_interaction,".\n",sep=""))  
               
        # Results for user defined model effects values
        if (!is.null(keepList)){
            # Model fit is selected according to user defined model effects    
            if(user_keep_batch && user_keep_equalvar){
                # Model 1
                model= model_MM
            }else if(user_keep_batch && !user_keep_equalvar){
                # Model 2
                model= model_hetvariance  
            }else if(!user_keep_batch && user_keep_equalvar){
                # Model 1A
                model= model_withoutbatch
            }else if(!user_keep_batch && !user_keep_equalvar){
                # Modify model 1A to heterogeneous residual variances
                model= do.call("gls", args=list(model.formula, weights=varIdent(form=~1|Genotype), x, na.action="na.omit"))
            }
            
            if(numberofgenders==2){
                if (equation=="withWeight"){   
                    interactionTest=anova(model, type="marginal")$"p-value"[5]           
                }
                else{
                    interactionTest=anova(model, type="marginal")$"p-value"[4]      
                } 
            }
            else {
                interactionTest = NA
            }
            
            compList <- (keepList==c(keep_batch,keep_equalvar,keep_weight,keep_gender,keep_interaction))
            
            if (length(compList[compList==FALSE])>0 && outputMessages)
                message("Warning:\nCalculated values differ from user defined values for model effects.\n")
            
            keep_weight <- user_keep_weight
            keep_gender <- user_keep_gender
            keep_interaction <- user_keep_interaction
            keep_batch <- user_keep_batch
            keep_equalvar <- user_keep_equalvar
            
        }
        
        if (outputMessages)
            message(paste("Information:\nEquation: '",equation,"'.\n",sep=""))   

        
        result <- new("PhenTestResult",list(model.output=model,depVariable=depVariable,equation=equation,method="MM", 
                        model.effect.batch=keep_batch,model.effect.variance=keep_equalvar,model.effect.interaction=keep_interaction,
                        model.output.interaction=interactionTest,model.effect.gender=keep_gender,model.effect.weight=keep_weight,
                        numberGenders=numberofgenders,pThreshold=pThreshold))
        
    }
    else if (method=="FE") {
        # Fisher Exact Test 
        if (outputMessages)
            message(paste("Information:\nMethod: Fisher Exact Test framework.\n",sep="")) 
        
        result <- buildFisherExactTest(phenList,depVariable,outputMessages)
    }
    else {
        if (outputMessages)   
            message(paste("Error:\nMethod define in the 'method' argument '",method,"' is not supported.\nAt the moment we are supporting 'MM' value for Mixed Model framework and 'FE' value for Fisher Exact Test framework.\n",sep=""))
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))      
        stop()
        
    }        
    
    return(result)   
}

#-----------------------------------------------------------------------------------
modelFormula <- function(equation, numberofgenders, depVariable)
# Create formula for the start model based on equation and number of genders in the data
{
    
    model.formula <- switch(equation,
            # Eq.2
            withWeight = {
                # Fixed effects: 1) Genotype 2) Gender 3) Genotype by Gender interaction 4) Weight
                if(numberofgenders==2){                            
                    as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Genotype*Gender","Weight", sep= "+")))
                }else{ 
                    as.formula(paste(depVariable, "~", paste("Genotype", "Weight",  sep= "+"))) 
                }         
            },
            # Eq.1 
            withoutWeight = {
                # Fixed effects: 1) Genotype 2) Gender 3) Genotype by Gender interaction
                if(numberofgenders==2){
                    as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Genotype*Gender", sep= "+")))
                }else{ 
                    as.formula(paste(depVariable, "~", paste("Genotype",  sep= "+"))) 
                } 
            }
            )
    return(model.formula)
}
