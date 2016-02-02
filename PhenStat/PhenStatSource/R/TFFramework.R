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
## TFFramework.R contains functions for regression analysis of data with no more than 5 batches (assays dates)
## Time as a fixed effect
##------------------------------------------------------------------------------
## Cleans dataset to make it suitable for the TF framework - 
## data points in all genotype/batch level combinations at least for one sex
TFDataset <- function(phenList=NULL, depVariable=NULL, outputMessages=TRUE, forDecisionTree=FALSE)
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
        # extract data frame 
        x <- getDataset(phenList)
        columnOfInterest <- getColumn(phenList,depVariable)
        # 3
        if (!batchIn(phenList)){
            stop_message <- paste("Error:\nBatch column is needed in the dataset to create make it suitable ",
                    "for the TF framework.\n",sep="")
        }
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
        # Remove NA records (Batch and/or depVariable are not defined)
        x<-x[!is.na(x$Batch),] 
        x<-x[!is.na(x[,c(depVariable)]),] 
        
        # Lists possible combinations
        Genotype_levels <- levels(factor(x$Genotype))
        Batch_levels <- levels(factor(x$Batch))
        Sex_levels <- levels(factor(x$Sex)) 
        #numberofSexes <- length(Sex_levels)
        tableWidth <- length(Sex_levels)*2+1
        countsAll <- nrow(x)
        countsRemovedRef <- 0
        countsRemovedTest <- 0
        cellWidth <- 13
        if (outputMessages){
            message("")
            if (outputMessages){
                message(paste("Data points containing '",depVariable,"' by batch levels:",sep=""))
            }
            line <- "-----------"
            message(printTabStyle(c(rep(line,each=tableWidth)),cellWidth))
            printGenotypes <- c("")
            for (j in 1:length(Genotype_levels)){   
                for (s in 1:length(Sex_levels)){
                    printGenotypes <- c(printGenotypes,Genotype_levels[j])
                }
            }
            message(printTabStyle(printGenotypes,cellWidth))
            message(printTabStyle(c(rep(line,each=tableWidth)),cellWidth))
            message(printTabStyle(c("Batch",Sex_levels,Sex_levels),cellWidth))
            message(printTabStyle(c(rep(line,each=tableWidth)),cellWidth))
        }
        removeRecords <- FALSE
        # Batch loop
        for (i in 1:length(Batch_levels)){
            BatchSubset <- subset(x, x$Batch==Batch_levels[i])
            sex_counts_batch <- c()
            removeBatch <- FALSE
            # Genotype loop
            for (j in 1:length(Genotype_levels)){           
                GenotypeBatchSubset <- subset(BatchSubset, 
                        BatchSubset$Genotype==Genotype_levels[j]) 
                # Sex loop
                for (s in 1:length(Sex_levels)){
                    GenotypeBatchSexSubset <- subset(GenotypeBatchSubset, 
                            GenotypeBatchSubset$Sex==Sex_levels[s]) 
                    columnOfInterestSubsetbySex <- na.omit(GenotypeBatchSexSubset[,c(depVariable)])
                    sex_counts_batch <- c(sex_counts_batch,length(columnOfInterestSubsetbySex) )
                }
                columnOfInterestSubset <- na.omit(GenotypeBatchSubset[,c(depVariable)])
                  
                if (length(columnOfInterestSubset)==0){
                    countsMinus <- nrow(x[x$Batch==Batch_levels[i],])
                    if (Genotype_levels[j]==testGenotype(phenList)){
                        countsRemovedRef <- countsRemovedRef + countsMinus 
                    }
                    else {
                        countsRemovedTest <- countsRemovedTest + countsMinus 
                    }
                    #x<- subset(x,!(x$Batch==Batch_levels[i]))
                    removeBatch <- TRUE
                    removeRecords <- TRUE
                }
            }   
            printBatchLevel <- Batch_levels[i]
            if (removeBatch){
                x<- subset(x,!(x$Batch==Batch_levels[i]))
                printBatchLevel <- paste("*", printBatchLevel)
            }
            if (outputMessages){
                message(printTabStyle(c(printBatchLevel,sex_counts_batch),cellWidth)) 
                message(printTabStyle(c(rep(line,each=tableWidth)),cellWidth))      
            }        
        }
        if (outputMessages && removeRecords){
            message("* - removed record(s)")
        }
        if (outputMessages || forDecisionTree){
            message("")
            if (forDecisionTree){
                message("TF Dataset Construction:")
            }
            message(paste("Number of batch levels left: ",length(levels(factor(x$Batch))),sep=""))
            message(paste("Records removed (reference genotype): ",round(countsRemovedRef*100/countsAll),"%",sep=""))
            message(paste("Records removed (test genotype): ",round(countsRemovedTest*100/countsAll),"%\n",sep=""))
        }
        

        if (length(levels(factor(x$Batch))) >= 2 && length(levels(factor(x$Batch))) <= 5) {
            new_phenList <- new("PhenList",datasetPL=x,
                    refGenotype = refGenotype(phenList),
                    testGenotype = testGenotype(phenList),
                    hemiGenotype = hemiGenotype(phenList))
            
            #PhenList(dataset=x, testGenotype=phenList$testGenotype, 
             #       refGenotype=phenList$refGenotype, hemiGenotype=phenList$hemiGenotype,outputMessages=FALSE)
            return (new_phenList)
        }
        else {
            stop_message <- paste("Error: \nThere are not enough records for TF method after",
                    " the removal of all non-concurrent batches.",sep="")
            
            if (outputMessages)  { 
                message(stop_message)
            }
            return(phenList)
        }
    
        
    }
}
##------------------------------------------------------------------------------
## Formats list of words into 'tab delimented' string 
printTabStyle <- function(textList,positions,tabSep="|"){
    result <- ""
    for (j in 1:length(textList)){
                textToPrint <- textList[j]
                if (nchar(textToPrint)<positions){
                    empty <- positions - nchar(textToPrint)
                }
                else {
                    empty <- 1
                }
                for (i in 1:empty){
                    textToPrint <- paste(" ",textToPrint,sep="")
                }
                result <- paste(result,tabSep,textToPrint," ",sep="")
    }
    return (paste(result,tabSep,sep=""))
}
##------------------------------------------------------------------------------
## Creates formula for the start model based on equation and number of Sexes
## in the data
startTFModel <- function(phenList, depVariable, equation="withWeight",
        outputMessages=TRUE, pThreshold=0.05, keepList=NULL)
{
    x <- getDataset(phenList)
    
    #if (!is.null(keepList)){
        ## User's values for effects
        #user_keep_weight <- keepList[3]
        #user_keep_sex <- keepList[4]
        #user_keep_interaction <- keepList[5]
        #user_keep_batch <- keepList[1]
        #user_keep_equalvar <- keepList[2]
        
        #if (!('Weight' %in% colnames(x)) && user_keep_weight){
        #    if (outputMessages)
        #    message(paste("Warning:\nWeight column is missed in the dataset.",
        #            "'keepWeight' is set to FALSE.\n"))
            
        #    user_keep_weight <- FALSE
        #}
        
        
        #if (equation=="withoutWeight" && user_keep_weight){
        #    if (outputMessages)
        #    message(paste("Warning:\nEquation is set to 'withoutWeight'.",
        #                    "'keepWeight' is set to FALSE.\n"))
            
        #    user_keep_weight <- FALSE
        #}
        
        #if (!('Batch' %in% colnames(x)) && user_keep_batch){
        #    if (outputMessages)
        #    message(paste("Warning:\nBatch column is missed in the dataset.",
        #            "'keepBatch' is set to FALSE.\n"))
            
        #    user_keep_batch <- FALSE
        #}
        
    #}
    
    numberofSexes <- length(levels(x$Sex))
    # Averages for percentage changes - is the ratio of the genotype effect for a sex relative to 
    # the wildtype signal for that variable for that sex - calculation        
    WT <- subset(x,x$Genotype==refGenotype(phenList))
    #mean_all <- mean(WT[,c(depVariable)],na.rm=TRUE)  
    mean_all <- mean(x[,c(depVariable)],na.rm=TRUE)  
    mean_list <- c(mean_all)  
    if (numberofSexes==2){  
      #WT_f <- subset(WT,WT$Sex=="Female")
      #WT_m <- subset(WT,WT$Sex=="Male")
      #mean_f <- mean(WT_f[,c(depVariable)],na.rm=TRUE)
      #mean_m <- mean(WT_m[,c(depVariable)],na.rm=TRUE)
      all_f <- subset(x,x$Sex=="Female")
      all_m <- subset(x,x$Sex=="Male")
      mean_f <- mean(all_f[,c(depVariable)],na.rm=TRUE)
      mean_m <- mean(all_m[,c(depVariable)],na.rm=TRUE)
      mean_list <- c(mean_all,mean_f,mean_m)  
    }
    # end of percentage change calculations 
    
    
    ## Start model formula: homogenous residual variance,
    ## genotype and sex interaction included
    
    model.formula  <- switch(equation,
                ## Eq.2
                withWeight = {
                    ## Fixed effects: 1) Genotype 2) Sex 3) Genotype by Sex
                    ## interaction 4) Weight
                    if(numberofSexes==2){
                        as.formula(paste(depVariable, "~", paste("Genotype",
                                                "Sex", "Genotype*Sex","Weight", sep= "+")))
                    }else{
                        as.formula(paste(depVariable, "~", paste("Genotype",
                                                "Weight", sep= "+")))
                    }
                },
                ## Eq.1
                withoutWeight = {
                    ## Fixed effects: 1) Genotype 2) Sex 3) Genotype by Sex
                    ## interaction
                    if(numberofSexes==2){
                        as.formula(paste(depVariable, "~",
                                        paste("Genotype", "Sex", "Genotype*Sex", sep= "+")))
                    }else{
                        as.formula(paste(depVariable, "~",
                                        paste("Genotype", sep= "+")))
                    }
                }
                )
    model.noBatch <- do.call("gls",args=list(model.formula, data = x, na.action="na.omit", method="ML"))
    
    if ('Batch' %in% colnames(x)){
        x<-x[!is.na(x$Batch),] 
        
        model.formula.withBatch  <- switch(equation,
                    ## Eq.2
                    withWeight = {
                        ## Fixed effects: 1) Genotype 2) Sex 3) Genotype by Sex
                        ## interaction 4) Weight
                        if(numberofSexes==2){
                            as.formula(paste(depVariable, "~", 
                                            paste("Genotype", "Sex", "Genotype*Sex", "Weight", "Batch", sep= "+")))
                        }else{
                            as.formula(paste(depVariable, "~", 
                                            paste("Genotype", "Weight", "Batch", sep= "+")))
                        }
                    },
                    ## Eq.1
                    withoutWeight = {
                        ## Fixed effects: 1) Genotype 2) Sex 3) Genotype by Sex
                        ## interaction
                        if(numberofSexes==2){
                            as.formula(paste(depVariable, "~",
                                            paste("Genotype", "Sex", "Genotype*Sex", "Batch", sep= "+")))
                        }else{
                            as.formula(paste(depVariable, "~",
                                            paste("Genotype", "Batch", sep= "+")))
                        }
                    }
            )
        
        
        model.withBatch <- do.call("gls",args=list(model.formula.withBatch, data = x, na.action="na.omit", method="ML"))
        ## Result of the test for Batch significance (fixed effect 5.)
        #anova_results <- anova(model.withBatch, type="marginal")$"p-value" < pThreshold
        p.value.batch <- (anova(model.withBatch, model.noBatch)$p[2])/2
        
        keep_batch <- p.value.batch < pThreshold
        
        #if(numberofSexes==2){
        #    keep_batch <- anova_results[5]
        #}
        #else {
        #    keep_batch <- anova_results[4]
        #}
        
        if (keep_batch){
            model <- model.withBatch
            model.formula <- model.formula.withBatch
        }
        else {
            model <- model.noBatch
        }
    }
    else {
        ## No Batch effects
        keep_batch <- FALSE
        
        model <- model.noBatch       
        
    }
    
    ## MM fit of model formula with heterogeneous residual variances for
    ## genotype groups
    ## Model 1 assumes homogeneous residual variances
    ## Model 2 with heterogeneous residual variances
    model_hetvariance <-
    tryCatch(
            model_hetvariance <- do.call("gls", args=list(model.formula, x,
                            weights=varIdent(form=~1|Genotype), na.action="na.omit", method="ML")),
            error=function(error_mes) {
                if (outputMessages)
                message(paste("Warning:\nMixed model with heterogeneous ",
                                "residual variances for genotype groups is not ",
                                "fitting - false convergence.\nMixed model with ",
                                "homogeneous residual variances is used instead.\n",sep=""))
                
                model_hetvariance <- NULL
            }
            )
    
    if (!is.null(model_hetvariance)) {
        ## Test: the variance of the residuals is the same (homogeneous)
        ## for all genotype groups
        ## Hypothesis 2
        ## Null Hypothesis: all residual variances are equal
        ## Alternative Hypothesis: the residue variance is not equal
        p.value.variance <- (anova(model, model_hetvariance)$p[2])
        ## The result of the test for Hypothesis 2 will help to select a
        ## covariance structure for the residuals
        keep_equalvar <- p.value.variance>pThreshold
    }
    else {
        keep_equalvar <- TRUE
        #if (!is.null(keepList)){
         #   if (outputMessages && !user_keep_equalvar)
         #   message(paste("Warning:\n'keepEqualVariance' is set to TRUE ", 
         #                   "otherwise the model can't be fitted - false convergence.\n",sep=""))
            
         #   user_keep_equalvar <- TRUE
        #}
    }
    
    
    
    ## Model fit is selected according to test results
    if(!keep_equalvar){
        model <- model_hetvariance
    }
    
    
    ## Tests for significance of fixed effects using TypeI F-test from anova
    ## functionality by using selected model
    anova_results <- anova(model, type="marginal")$"p-value" < pThreshold
    
    #print(anova_results)
    #print(anova(model, type="marginal"))
    
    if(numberofSexes==2){
        ## Result of the test for Sex significance (fixed effect 1.)
        keep_sex <- anova_results[3]
        ## Eq.2
        if (equation=="withWeight"){
            ## Result of the test for weight significance  (fixed effect 3.)
            
            keep_weight <- anova_results[4]
            ## Result of the test for genotype by Sex interaction
            ## significance (fixed effect 2.)
            if (keep_batch) {
                keep_interaction <- anova_results[6]
                ## Technical results needed for the output
                ## Interaction test results are kept for the output
                interactionTest <- anova(model, type="marginal")$"p-value"[6]
            }
            else {
                keep_interaction <- anova_results[5]
                ## Technical results needed for the output
                ## Interaction test results are kept for the output
                interactionTest <- anova(model, type="marginal")$"p-value"[5]
            }
            
            
            
            
            
        }
        ## Eq.1
        else{
            ## Result of the test for weight significance  (fixed effect 3.)
            ## It's FALSE since here the equation 1 is used - without weight
            ## effect
            keep_weight <- FALSE
            ## Result of the test for genotype by Sex interaction
            ## significance (fixed effect 2.)
            if (keep_batch) {                
                keep_interaction <- anova_results[5]
                ## Interaction test results are kept for the output
                interactionTest <- anova(model, type="marginal")$"p-value"[5]
            }
            else {                
                keep_interaction <- anova_results[4]
                ## Interaction test results are kept for the output
                interactionTest <- anova(model, type="marginal")$"p-value"[4]
            }
            
        }
    }
    else {
        keep_sex <- FALSE
        keep_interaction <- FALSE
        interactionTest <- NA
        if (equation=="withWeight")
            keep_weight <- anova_results[3]
        else
            keep_weight <- FALSE
    }
    
    if (!keep_weight && equation=="withWeight") {
        equation="withoutWeight"
        if (outputMessages)
        message(paste("Since weight effect is not significant the equation ",
                        "'withoutWeight' should be used instead.",sep=""))
    }
    
    if (outputMessages)
    message(paste("Information:\nEquation: '",equation,"'.\n",sep=""))
    
    if (outputMessages)
    message(paste("Information:\nCalculated values for model effects are: ",
                    "keepBatch=",keep_batch,
                    ", keepEqualVariance=",keep_equalvar,
                    ", keepWeight=",keep_weight,
                    ", keepSex=",keep_sex,
                    ", keepInteraction=",keep_interaction,".\n",sep=""))
    
    ## Results for user defined model effects values
    #if (!is.null(keepList)){
    #    if (outputMessages)
    #    message(paste("Information:\nUser's values for model effects are: ",
    #                    "keepBatch=",user_keep_batch,
    #                    ", keepEqualVariance=",user_keep_equalvar,
    #                    ", keepWeight=",user_keep_weight,
    #                    ", keepSex=",user_keep_sex,
    #                    ", keepInteraction=",user_keep_interaction,".\n",sep=""))
        ## Model fit is selected according to user defined model effects
    #    if(user_keep_batch && user_keep_equalvar){
            ## Model 1
     #       model <- model.withBatch
     #   }else if(user_keep_batch && !user_keep_equalvar){
            ## Model 2
     #       model <- model_hetvariance
     #   }else if(!user_keep_batch && user_keep_equalvar){
            ## Model 1A
     #       model <- model.noBatch
     #   }else if(!user_keep_batch && !user_keep_equalvar){
            ## Modify model 1A to heterogeneous residual variances
     #       model <- do.call("gls", args=list(modelFormulaNoBatch(equation,numberofSexes, depVariable),
     #                      weights=varIdent(form=~1|Genotype), x, na.action="na.omit"))
     #   }
        
      #  if(numberofSexes==2){
      #      index <- 6
      #      if (equation=="withoutWeight"){
      #         index <- index-1
      #          
      #      }
      #      if (!keep_batch) {
      #          index <- index-1
      #      }
      #      interactionTest <- anova(model, type="marginal")$"p-value"[index]
      #  }
      #  else {
      #      interactionTest <- NA
      #  }
        
      #  compList <- (keepList==c(keep_batch,keep_equalvar,keep_weight,
      #                  keep_sex,keep_interaction))
        
      #  if (length(compList[compList==FALSE])>0 && outputMessages)
      #  message("Warning:\nCalculated values differ from user defined values for model effects.\n")
        
      #  keep_weight <- user_keep_weight
      #  keep_sex <- user_keep_sex
      #  keep_interaction <- user_keep_interaction
      #  keep_batch <- user_keep_batch
      #  keep_equalvar <- user_keep_equalvar
        
    #}
    
    # TF modelling output
    linearRegressionOutput <- list(model.output=model,
        equation=equation,
        model.effect.batch=keep_batch,
        model.effect.variance=keep_equalvar,
        model.effect.interaction=keep_interaction,
        model.output.interaction=interactionTest,
        model.effect.sex=keep_sex,
        model.effect.weight=keep_weight,
        numberSexes=numberofSexes,
        pThreshold=pThreshold,
        model.formula.genotype=model.formula,
        model.output.averageRefGenotype=mean_list)
        

    result <- new("PhenTestResult",
            analysedDataset=x,
            depVariable=depVariable,
            refGenotype=refGenotype(phenList),
            testGenotype=testGenotype(phenList),
            method="TF",
            transformationRequired = FALSE,
            lambdaValue = integer(0),
            scaleShift = integer(0),
            analysisResults=linearRegressionOutput)

    return(result)
}

##------------------------------------------------------------------------------
## Works with PhenTestResult object created by testDataset function.
## Builds final model based on the significance of different model effects,
## depVariable and equation stored in phenTestResult object (see testDataset.R).
finalTFModel <- function(phenTestResult, outputMessages=TRUE)
{
    ## Checks and stop messages
    stop_message <- ""
    
    ## Check PhenTestResult object
    if(is(phenTestResult,"PhenTestResult")) {

        finalResult <- phenTestResult
        linearRegressionOutput <- analysisResults(phenTestResult)
        x <- analysedDataset(finalResult)
        depVariable <- getVariable(finalResult)
        equation <- linearRegressionOutput$equation
        keep_weight <- linearRegressionOutput$model.effect.weight
        keep_sex <- linearRegressionOutput$model.effect.sex
        keep_interaction <- linearRegressionOutput$model.effect.interaction
        keep_batch <- linearRegressionOutput$model.effect.batch
        keep_equalvar <- linearRegressionOutput$model.effect.variance

        ## Stop function if there are no datasets to work with
        if(is.null(x)){
            stop_message <- "Error:\nPlease create a PhenList object first and run function 'testDataset'.\n"
        }
        
        ## Stop function if there are no enough input parameters
        if (is.null(equation) || is.null(depVariable) || is.null(keep_batch) || is.null(keep_equalvar)
            || is.null(keep_sex) || is.null(keep_interaction)){
            stop_message <- "Error:\nPlease run function 'testDataset' first.\n"
        }
    }else{
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
    
    numberofSexes <- linearRegressionOutput$numberSexes
    
    
    ## Build final null model
    ## Goal:  to test fixed effects of the model and based on the output
    ## build the final null model formula for later
    ## testing - as a null model it automatically excludes genotype and
    ## interaction term.
    ## The anova function tests the fixed effects associated by treatment
    ## with a null hypothesis that the regression
    ## coefficients are equal to zero  and an alternative hypothesis that the
    ## regression coefficient are not equal to zero.
    ## If the p-values of these tests are less than 0.05 we reject the null
    ## and accept the alternative that the are
    ## significant components of the model and should be included.
    ## If no terms are significant a model can be build with just an
    ## intercept element this is specified as
    ## "model.formula <- as.formula(paste(depVariable, "~", "1"))"
    
    ## Null model: genotype is not significant
    batch_formula <- ""
    if (keep_batch){
        batch_formula <- "+ Batch"
    }
    model_null.formula <- switch(equation,
        withWeight = {
            ## Eq.2
            if(numberofSexes==2){
                if(!keep_sex){
                    as.formula(paste(depVariable, "~", "Weight",batch_formula))
                }else{
                    as.formula(paste(depVariable, "~", 
                    paste("Sex", "Weight", sep= "+"),batch_formula))
                }
            }else{
                as.formula(paste(depVariable, "~", "Weight",batch_formula))
            }
        },
        withoutWeight = {
            ## Eq.1
            if(numberofSexes==2){
                if(!keep_sex && !keep_interaction){
                    if (keep_batch){
                        as.formula(paste(depVariable, "~", "Batch"))
                    }
                    else {
                        as.formula(paste(depVariable, "~", "1"))
                    }
                }else{
                    as.formula(paste(depVariable, "~", "Sex",batch_formula))
                }
            }else{
                if (keep_batch){
                    as.formula(paste(depVariable, "~", "Batch"))
                }
                else {
                    as.formula(paste(depVariable, "~", "1"))
                }
            }
        }
    )
    

    ## Alternative model: genotype is significant
    model_genotype.formula <- switch(equation,
        withWeight = {
            ## Eq.2
            if(numberofSexes==2){
                if ((keep_sex && keep_weight && keep_interaction)|
                        (!keep_sex && keep_weight && keep_interaction)){
                    as.formula(paste(depVariable, "~",
                    paste("Sex", "Genotype:Sex", "Weight", sep= "+"),batch_formula))
                
                } else if(keep_sex && keep_weight && !keep_interaction){
                    as.formula(paste(depVariable, "~",
                    paste("Genotype", "Sex", "Weight", sep= "+"),batch_formula))
                
                } else if(!keep_sex && keep_weight && !keep_interaction){
                    as.formula(paste(depVariable, "~",
                    paste("Genotype", "Weight", sep= "+"),batch_formula))
                }
            
            }else{
                as.formula(paste(depVariable, "~",
                paste("Genotype", "Weight", sep="+"),batch_formula))
            }
        },
        withoutWeight = {
        ## Eq.1
            if(numberofSexes==2){
                if (!keep_sex  && !keep_interaction){
                    as.formula(paste(depVariable, "~", "Genotype",batch_formula))
                } else if((keep_sex && keep_interaction)|(!keep_sex && keep_interaction)){
                    as.formula(paste(depVariable, "~",
                    paste("Sex",  "Genotype:Sex", sep= "+"),batch_formula))
                } else if(keep_sex && !keep_interaction){
                    as.formula(paste(depVariable, "~",
                    paste("Genotype","Sex", sep= "+"),batch_formula))
                }
            }else{
                as.formula(paste(depVariable, "~", "Genotype",batch_formula))
            }
        }
    )
    
    #print(model_genotype.formula)
    #print(model_null.formula)
    ## Test: genotype groups association with dependent variable
    ## Null Hypothesis: genotypes are not associated with dependent variable
    ## Alternative Hypothesis: genotypes are associated with dependent
    ## variable
    if(keep_equalvar){
        model_genotype <- do.call("gls", args = list(model_genotype.formula,
        x, method="ML", na.action="na.omit"))
        
        model_null <- do.call("gls", args=list(model_null.formula, x,
        method="ML", na.action="na.omit"))
        
        p.value <- (anova(model_genotype, model_null)$p[2])
    }else {
        model_genotype <- do.call("gls", args = list(model_genotype.formula,
        x,weights=varIdent(form=~1|Genotype),method="ML", na.action="na.omit"))
        
        model_null <- do.call("gls", args=list(model_null.formula,
        x,weights=varIdent(form=~1|Genotype),method="ML", na.action="na.omit"))
        
        p.value <- (anova(model_genotype, model_null)$p[2])
    }
    
    ## Final model version with na.exclude and REML method
    if(keep_equalvar){
        model_genotype <- do.call("gls", args = list(model_genotype.formula,
        x, na.action="na.exclude"))
    
    }else {
        model_genotype <- do.call("gls", args = list(model_genotype.formula,
        x,weights=varIdent(form=~1|Genotype), na.action="na.exclude"))
    }
    

    ## Store the results after the final TF modelling step
    linearRegressionOutput$model.output <- model_genotype
    linearRegressionOutput$model.null <- model_null
    linearRegressionOutput$model.output.genotype.nulltest.pVal <- p.value
    linearRegressionOutput$model.formula.null <- model_null.formula
    linearRegressionOutput$model.formula.genotype <- model_genotype.formula
    #linearRegressionOutput$model.effect.variance <- keep_equalvar
    
    
    
    ## Parse modeloutput and choose output depending on model
    linearRegressionOutput$model.output.summary <- parserOutputTFSummary(linearRegressionOutput)

    # Percentage changes - is the ratio of the genotype effect for a sex relative to 
    # the wildtype signal for that variable for that sex - calculation   
    mean_list <- linearRegressionOutput$model.output.averageRefGenotype
    denominator <- mean_list[1]
    if(numberofSexes==2){
        # without weight
        if (is.na(linearRegressionOutput$model.output.summary['weight_estimate'])){ 
            if (!is.na(linearRegressionOutput$model.output.summary['sex_estimate']) &&
            !is.na(linearRegressionOutput$model.output.summary['sex_FvKO_estimate']))
            {
                denominator_f <- linearRegressionOutput$model.output.summary['intercept_estimate']
                denominator_m <- linearRegressionOutput$model.output.summary['intercept_estimate']+
                linearRegressionOutput$model.output.summary['sex_estimate']
                denominator_f <- denominator
                denominator_m <- denominator
                ratio_f <- linearRegressionOutput$model.output.summary['sex_FvKO_estimate']/denominator_f                       
                ratio_m <- linearRegressionOutput$model.output.summary['sex_MvKO_estimate']/denominator_m 
            }
            else if (!is.na(linearRegressionOutput$model.output.summary['sex_FvKO_estimate']))
            {
                #denominator <- linearRegressionOutput$model.output.summary['intercept_estimate']
                ratio_f <- linearRegressionOutput$model.output.summary['sex_FvKO_estimate']/denominator                            
                ratio_m <- linearRegressionOutput$model.output.summary['sex_MvKO_estimate']/denominator            
            }
            else if (!is.na(linearRegressionOutput$model.output.summary['sex_estimate']))
            {
                denominator_f <- linearRegressionOutput$model.output.summary['intercept_estimate']
                denominator_m <- linearRegressionOutput$model.output.summary['intercept_estimate']+
                linearRegressionOutput$model.output.summary['sex_estimate']
                denominator_f <- denominator
                denominator_m <- denominator
                ratio_f <- linearRegressionOutput$model.output.summary['genotype_estimate']/denominator_f                        
                ratio_m <- linearRegressionOutput$model.output.summary['genotype_estimate']/denominator_m 
            }
            else
            {
                #denominator <- linearRegressionOutput$model.output.summary['intercept_estimate']
                ratio_f <- linearRegressionOutput$model.output.summary['genotype_estimate']/denominator                            
                ratio_m <- ratio_f                      
            }
        }
        # with weight
        else{
            mean_list <- linearRegressionOutput$model.output.averageRefGenotype
            denominator_f <- mean_list[2]
            denominator_m <- mean_list[3]
            denominator_f <- denominator
            denominator_m <- denominator
            if (!is.na(linearRegressionOutput$model.output.summary['sex_estimate']) &&
            !is.na(linearRegressionOutput$model.output.summary['sex_FvKO_estimate']))
            {
                ratio_f <- linearRegressionOutput$model.output.summary['sex_FvKO_estimate']/denominator_f                       
                ratio_m <- linearRegressionOutput$model.output.summary['sex_MvKO_estimate']/denominator_m 
            }
            else if (!is.na(linearRegressionOutput$model.output.summary['sex_FvKO_estimate']))
            {
                ratio_f <- linearRegressionOutput$model.output.summary['sex_FvKO_estimate']/denominator_f                            
                ratio_m <- linearRegressionOutput$model.output.summary['sex_MvKO_estimate']/denominator_m            
            }
            else 
            {
                ratio_f <- linearRegressionOutput$model.output.summary['genotype_estimate']/denominator_f                        
                ratio_m <- linearRegressionOutput$model.output.summary['genotype_estimate']/denominator_m 
            }
        }
        linearRegressionOutput$model.output.percentageChanges <- c(ratio_f*100,ratio_m*100)
        names(linearRegressionOutput$model.output.percentageChanges) <- c('female*genotype ratio','male*genotype ratio')
    }
    else{
        # without weight
        if (is.na(linearRegressionOutput$model.output.summary['weight_estimate'])){ 
            #denominator <- linearRegressionOutput$model.output.summary['intercept_estimate']
            ratio_f <- linearRegressionOutput$model.output.summary['genotype_estimate']/denominator                            
        }
        # with weight
        else{
            mean_list <- linearRegressionOutput$model.output.averageRefGenotype
            denominator <- mean_list[1]
            ratio_f <- linearRegressionOutput$model.output.summary['genotype_estimate']/denominator   
        }
        
    linearRegressionOutput$model.output.percentageChanges <- c(ratio_f*100)
    names(linearRegressionOutput$model.output.percentageChanges) <- c('all*genotype ratio')
    }
    # end of percentage changes calculation
    finalResult@analysisResults <- linearRegressionOutput 
    ## Assign MM quality of fit
    finalResult@analysisResults$model.output.quality <- testFinalModel(finalResult)
    
   
    
    return(finalResult)
}
##------------------------------------------------------------------------------
## Parser model output summary and return in readable vector format
parserOutputTFSummary<-function(linearRegressionOutput)
{
    result <- linearRegressionOutput
    modeloutput_summary <- summary(result$model.output)
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

    lengthoftable <- {
        table_length <- NA
        
        if (result$equation=="withWeight"){
            if(result$numberSexes==2){
                if((result$model.effect.sex
                                && result$model.effect.interaction)|
                        (!result$model.effect.sex && result$model.effect.interaction)){
                    table_length <- 5
                }else if(result$model.effect.sex &&
                        !result$model.effect.interaction){
                    table_length <- 4
                }else {
                    table_length <- 3
                }
            }else{
                table_length <- 3
            }
        }
        ## Eq.1
        else{
            if(result$numberSexes==2){
                if((result$model.effect.sex
                                && result$model.effect.interaction)|
                        (!result$model.effect.sex && result$model.effect.interaction)){
                    table_length <- 4
                }else if(!result$model.effect.sex &&
                        !result$model.effect.interaction){
                    table_length <- 2
                }else{
                    table_length <- 3
                }
            }else{
                table_length <- 2
            }
        }
        
        table_length
    }
    
    if (result$model.effect.batch){
        lengthofbatch <- length(grep("Batch",row.names(modeloutput_summary[["tTable"]])))
        #lengthofbatch <- length(levels(result$model.dataset$Batch)) - 1
        #if (lengthofbatch == -1){
         #   lengthofbatch <- 1
        #}
        lengthoftable <- lengthoftable + lengthofbatch
    }
    # 1) intercept
    intercept_estimate = modeloutput_summary[["tTable"]][[1]]
    intercept_estimate_SE = modeloutput_summary[["tTable"]][[(1+lengthoftable)]]

    if((result$model.effect.sex && result$model.effect.interaction)
    |( !result$model.effect.sex && result$model.effect.interaction)){
        sex_index <- match(c("SexMale"),row.names(modeloutput_summary[["tTable"]]))
        #sex_FvKO_index <- lengthoftable-1
        #sex_MvKO_index <- lengthoftable
        sex_FvKO_index <- grep("SexFemale:Genotype",row.names(modeloutput_summary[["tTable"]]))[1]
        sex_MvKO_index <- grep("SexMale:Genotype",row.names(modeloutput_summary[["tTable"]]))[1]
        if (is.na(sex_index)){
            sex_index <- match(c("SexFemale"),row.names(modeloutput_summary[["tTable"]]))
            #sex_FvKO_index <- lengthoftable
            #sex_MvKO_index <- lengthoftable-1
        }

        sex_estimate <- modeloutput_summary[["tTable"]][[sex_index]]
        sex_estimate_SE <- modeloutput_summary[["tTable"]][[(sex_index+lengthoftable)]]
        sex_p_value <- modeloutput_summary[["tTable"]][[(sex_index+3*lengthoftable)]]
        
        sex_FvKO_estimate= modeloutput_summary[["tTable"]][[sex_FvKO_index]]
        sex_FvKO_SE=modeloutput_summary[["tTable"]][[(sex_FvKO_index+lengthoftable)]]
        sex_FvKO_p_value=modeloutput_summary[["tTable"]][[(sex_FvKO_index+3*lengthoftable)]]
        sex_MvKO_estimate=modeloutput_summary[["tTable"]][[sex_MvKO_index]]
        sex_MvKO_SE=modeloutput_summary[["tTable"]][[(sex_MvKO_index+lengthoftable)]]
        sex_MvKO_p_value=modeloutput_summary[["tTable"]][[(sex_MvKO_index+3*lengthoftable)]]   
    
    } else if( !result$model.effect.sex && !result$model.effect.interaction){    

        genotype_index <- grep("Genotype",row.names(modeloutput_summary[["tTable"]]))[1] 
        genotype_estimate = modeloutput_summary[["tTable"]][[genotype_index]]
        genotype_estimate_SE = modeloutput_summary[["tTable"]][[(genotype_index+lengthoftable)]]
        genotype_p_value =  modeloutput_summary[["tTable"]][[(genotype_index+3*lengthoftable)]]
    
    }else{
        sex_index <- match(c("SexMale"),row.names(modeloutput_summary[["tTable"]]))
        if (is.na(sex_index)){
            sex_index <- match(c("SexFemale"),row.names(modeloutput_summary[["tTable"]]))
        }
        genotype_index <- grep("Genotype",row.names(modeloutput_summary[["tTable"]]))[1]      
        genotype_estimate = modeloutput_summary[["tTable"]][[genotype_index]]
        genotype_estimate_SE = modeloutput_summary[["tTable"]][[(genotype_index+lengthoftable)]]
        genotype_p_value =  modeloutput_summary[["tTable"]][[(genotype_index+3*lengthoftable)]]
        sex_estimate=modeloutput_summary[["tTable"]][[sex_index]]
        sex_estimate_SE=modeloutput_summary[["tTable"]][[(sex_index+lengthoftable)]]
        sex_p_value= modeloutput_summary[["tTable"]][[(sex_index+3*lengthoftable)]]
    }

    if (result$equation=="withWeight"){
        weight_index <- match(c("Weight"),row.names(modeloutput_summary[["tTable"]]))
        weight_estimate=modeloutput_summary[["tTable"]][[weight_index]]
        weight_estimate_SE=modeloutput_summary[["tTable"]][[(weight_index+lengthoftable)]]
        weight_p_value=modeloutput_summary[["tTable"]][[(weight_index+3*lengthoftable)]]
    }

    output <- c(genotype_estimate, genotype_estimate_SE,  genotype_p_value,
            sex_estimate, sex_estimate_SE,  sex_p_value, 
            weight_estimate, weight_estimate_SE, weight_p_value, 
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