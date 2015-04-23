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
## simpleOutput.R contains simpleOutput, generateGraphs functions
##------------------------------------------------------------------------------
## Wrapper to prepare the output of the modelling and testing results in simple 
## user friendly form. Assumes that modelling results are stored in the
## phenTestResult object (output from functions testDataset and buildFinalModel)
summaryOutput <- function(phenTestResult,phenotypeThreshold=0.01)
{
    line <- "----------------------------------------------------------------------------"
    message(paste("\nTest for dependent variable:\n*** ",
                    getVariable(phenTestResult), transformationText(phenTestResult),
                    " ***",sep=""))
    
    message(paste("\nMethod:\n*** ",switch(method(phenTestResult),
                            MM = "Mixed Model framework",
                            FE = "Fisher Exact Test framework",
                            RR = "Reference Ranges Plus framework",
                            LR = "Logistic Regression",
                            TF = "Time as Fixed Effect framework")," ***",
                    "\n",sep=""))
    
    if (method(phenTestResult) %in% c("MM","TF","LR")) {
        linearRegressionOutput <- analysisResults(phenTestResult)
        
        message(line)
        message("Model Output")
        message(line)        
        message(paste("Final fitted model:",
                        format(linearRegressionOutput$model.formula.genotype)))
        message(paste("Was batch significant?",linearRegressionOutput$model.effect.batch))
        if (method(phenTestResult)!="LR")
            message(paste("Was variance equal?",linearRegressionOutput$model.effect.variance))
   
        printLROutput(phenTestResult,phenotypeThreshold)
        
        message(paste("\n",line,sep=""))
        message("Classification Tag")
        message(line)
        message(classificationTag(phenTestResult,phenotypeThreshold=phenotypeThreshold))
        
        message(paste("\n",line,sep=""))
        message("Model Output Summary")
        message(line)
        if (method(phenTestResult) %in% c("MM","TF")) 
            summary(linearRegressionOutput$model.output)$tTable
        else {
            error_estimates=sqrt(diag(vcov(linearRegressionOutput$model.output)))    
            ab=cbind(linearRegressionOutput$model.output$coefficients,error_estimates, 
                    linearRegressionOutput$model.output$ci.lower, 
                    linearRegressionOutput$model.output$ci.upper, 
                    linearRegressionOutput$model.output$prob)
            ab=as.data.frame(ab)
            colnames(ab)=c("Value", "Std.Error", "ci.lower", "ci.upper", "p-value")
            print(ab)
        }
    }
    
    else if (method(phenTestResult) %in% c("FE")){
        colnum <- 1
        message(line)
        message(paste("Model Output ('*' highlights results with p-values less than threshold ",phenotypeThreshold
                        ,")",sep=""))
        message(line)
        for (i in seq_along(analysisResults(phenTestResult))) {
            val <- analysisResults(phenTestResult)[[i]]
            message(cat(subsetText(val)))
            print(getMatrix(val),quote=FALSE)
            message("\n")
        }    
           
        
        message(paste("\n",line,sep=""))
        message("Classification Tag")
        message(line)
        
        message(classificationTag(phenTestResult,phenotypeThreshold=phenotypeThreshold))

        message(paste("\n",line,sep=""))
        message("Count Matrices")
        message(line)
        for (i in seq_along(analysisResults(phenTestResult))) {
            val <- analysisResults(phenTestResult)[[i]]
            message(cat(subsetText(val)))
            print(matrixCount(val))
            message("\n")
        } 
    }    
    
    else if (method(phenTestResult) %in% c("RR")){
        x <- analysedDataset(phenTestResult)
        noSexes <- length(levels(x$Sex))
        message(line)
        message("Model Output")
        message(line)
        colnum <- 1
        #cat("\n1) High vs Normal/Low\n")
        nl <- data.frame(nr=c(1,2,3))
        for (i in seq_along(analysisResults(phenTestResult))) {
            val <- analysisResults(phenTestResult)[[i]]
            if (comparison(val)=="High vs Normal/Low"){
                nl[paste("new",i)]<-as.vector(getColumnView(val))                  
            }
        }
        nl<-nl[ , -which(names(nl) %in% c("nr"))]
        nh <- data.frame(nr=c(1,2,3))
        for (i in seq_along(analysisResults(phenTestResult))) {
            val <- analysisResults(phenTestResult)[[i]]
            if (comparison(val)=="Low vs Normal/High"){
                nh[paste("new",i)]<-as.vector(getColumnView(val))                  
            }
        }
        nh<-nh[ , -which(names(nh) %in% c("nr"))]
        
        if (noSexes==2){
            colnames(nl) <- head(nl,1)
            nl <- nl[-1,]
            rownames(nl) <- c("p-value","ES")
            cat("\n1) High vs Normal/Low\n")
            print(nl,quote=FALSE)
            colnames(nh) <- head(nh,1)
            nh <- nh[-1,]
            rownames(nh) <- c("p-value","ES")
            cat("\n2) Low vs Normal/High\n")
            print(nh,quote=FALSE)
        }
        else {
            all <- cbind(nl,nh)
            all <- all[-1,]
            rownames(all) <- c("p-value","ES")
            colnames(all) <- c("High vs Normal/Low","Low vs Normal/High")
            print(all,quote=FALSE)
        }
        
        
        
        message(paste("\n",line,sep=""))
        message("Classification Tag")
        message(line)
        message(classificationTag(phenTestResult,phenotypeThreshold=phenotypeThreshold))
        
        message(paste("\n",line,sep=""))
        message("Thresholds")
        message(line)
        print(parameters(phenTestResult),quote=FALSE)
        message(paste("\n",line,sep=""))
        message("Count Matrices")
        message(line)
        for (i in seq_along(analysisResults(phenTestResult))) {
            val <- analysisResults(phenTestResult)[[i]]
            message(cat(subsetText(val)))
            print(matrixCount(val))
            message("\n")
        } 
        

    }  
}


##------------------------------------------------------------------------------
## Function for linear regression output
printLROutput <- function(phenTestResult,phenotypeThreshold=0.01)
{
    
    linearRegressionOutput <- analysisResults(phenTestResult)
    
    if (linearRegressionOutput$model.effect.interaction==TRUE){
        effectValuesMales <- c(linearRegressionOutput$model.output.summary["sex_MvKO_estimate"],
                linearRegressionOutput$model.output.summary["sex_MvKO_SE"])
        effectValuesFemales <- c(linearRegressionOutput$model.output.summary["sex_FvKO_estimate"],
                linearRegressionOutput$model.output.summary["sex_FvKO_SE"])
        
        if (phenTestResult@transformationRequired) {
            effectValuesMales <- as.numeric(reverseTransformValues(effectValuesMales,phenTestResult@lambdaValue,
                            phenTestResult@scaleShift))
            effectValuesFemales <- as.numeric(reverseTransformValues(effectValuesFemales,phenTestResult@lambdaValue,
                            phenTestResult@scaleShift))
        }
        effectValues <- c(effectValuesMales,effectValuesFemales)
    }
    else { 
        effectValues <- getGenotypeEffect(phenTestResult)  
    }
    
    message(paste("Genotype p-value:",
                    sprintf("%0.6e",linearRegressionOutput$model.output.genotype.nulltest.pVal)))
    
    original_scale <- ""
    if (phenTestResult@transformationRequired){    
        original_scale <- " (original scale)"
    }
    
    if (linearRegressionOutput$model.effect.interaction==FALSE){
        message(paste("Genotype effect",original_scale,": ",
                        sprintf("%.4f",round(effectValues[1],digits=4))," +/- ",
                        sprintf("%.4f",abs(round(effectValues[2],digits=4))),sep=""))
    }
    else {      
        message(paste("Genotype by male effect",original_scale,": ",
                        sprintf("%.4f",round(as.numeric(effectValues[1]),digits=4))," +/- ",
                        sprintf("%.4f",abs(round(as.numeric(effectValues[2]),digits=4))),sep=""))
        message(paste("Genotype by female effect",original_scale,": ",
                        sprintf("%.4f",round(as.numeric(effectValues[3]),digits=4))," +/- ",
                        sprintf("%.4f",abs(round(as.numeric(effectValues[4]),digits=4))),sep=""))

    }
    
    if (linearRegressionOutput$model.effect.interaction)
    {
        sexualDimorphism <- "yes"
    }
    else {
        sexualDimorphism <- "no"
    }
    message(paste("Was there evidence of sexual dimorphism? ",
                    sexualDimorphism," (p-value ",
                    sprintf("%0.6e",linearRegressionOutput$model.output.interaction),")",sep=""))
    
    if (!is.null(linearRegressionOutput$model.output.percentageChanges)){
        if (linearRegressionOutput$numberSexes==2){
            message(paste("Genotype percentage change Female:",
                            round(linearRegressionOutput$model.output.percentageChanges[1],digits=2)),"%")
            message(paste("Genotype percentage change Male:",
                            round(linearRegressionOutput$model.output.percentageChanges[2],digits=2)),"%")
        } 
        else {
            if ("Female" %in% levels(linearRegressionOutput$model.dataset$Sex)){
                message(paste("Genotype percentage change Female:",
                                round(linearRegressionOutput$model.output.percentageChanges[1],digits=2)),"%")
            }  
            else {
                message(paste("Genotype percentage change Male:",
                                round(linearRegressionOutput$model.output.percentageChanges[1],digits=2)),"%")
                
            }
        }   
    }
    
    

}   

##------------------------------------------------------------------------------
## Generate graphs and store them in the defined directory
generateGraphs <- function(phenTestResult, dir, graphingName=NULL, type="Xlib")
{
    if (method(phenTestResult) %in% c("MM","TF","LR")){
        
        if (is.null(graphingName)){
            graphingName = getVariable(phenTestResult)
        }
        
        
        #1
        graph_name=file.path(dir, "qqplotGenotype.png")
        png(graph_name,type=type)
        qqplotGenotype(phenTestResult,outputMessages=FALSE)
        dev.off()
        
        #2
        graph_name=file.path(dir, "plotResidualPredicted.png")
        png(graph_name,type=type)
        plotResidualPredicted(phenTestResult)
        dev.off()
        
        #PhenList object
        phenList <- new("PhenList",datasetPL=analysedDataset(phenTestResult),
                refGenotype = refGenotype(phenTestResult),
                testGenotype = testGenotype(phenTestResult))
        
        if (batchIn(phenTestResult)){
            #3
            graph_name=file.path(dir, "scatterplotSexGenotypeBatch.png")
            png(graph_name,type=type)
            scatterplotSexGenotypeBatch(phenList,
                    getVariable(phenTestResult), graphingName=graphingName)            
            dev.off()
            #4
            graph_name=file.path(dir, "boxplotResidualBatch.png")
            png(graph_name,type=type)
            boxplotResidualBatch(phenTestResult,outputMessages=FALSE)            
            dev.off()
        }
        
        if (weightIn(phenTestResult)){
            #5
            graph_name=file.path(dir, "scatterplotGenotypeWeight.png")
            png(graph_name,type=type)
            scatterplotGenotypeWeight(phenList,
                    getVariable(phenTestResult), graphingName=graphingName)
            
            dev.off()
        }
        
        # Two graphs for MM and Batch effect only
        if (batchIn(phenTestResult) && (method(phenTestResult)=="MM")){
            if (analysisResults(phenTestResult)$model.effect.batch){
                #6
                graph_name=file.path(dir, "qqplotRandomEffects.png")
                png(graph_name,type=type)
                qqplotRandomEffects(phenTestResult,outputMessages=FALSE)
                dev.off()
                #7
                graph_name=file.path(dir, "qqplotRotatedResiduals.png")
                png(graph_name,type=type)
                qqplotRotatedResiduals(phenTestResult,outputMessages=FALSE)
                dev.off()
            }
        }    
        
    }
    else {
        graph_name=file.path(dir, "categoricalBarPlot.png")
        png(graph_name,type=type)
        categoricalBarplot(phenTestResult,outputMessages=FALSE)        
        dev.off()
    } 
    
}    
##------------------------------------------------------------------------------