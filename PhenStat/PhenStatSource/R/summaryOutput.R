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
    message(paste("\nTest for dependent variable:\n*** ",phenTestResult$depVariable," ***",sep=""))
    
    message(paste("\nMethod:\n*** ",switch(phenTestResult$method,
                            MM = "Mixed Model framework",
                            FE = "Fisher Exact Test framework",
                            RR = "Reference Ranges Plus framework",
                            TF = "Time as Fixed Effect framework")," ***",
                    "\n",sep=""))
    
    if (phenTestResult$method %in% c("MM","TF")) {
        message(line)
        message("Model Output")
        message(line)
        
        message(paste("Genotype effect:",
                        sprintf("%.4f",round(phenTestResult$model.output.genotype.nulltest.pVal,digits=4))))
        message(paste("\nFinal fitted model:",
                        format(phenTestResult$model.formula.genotype)))
        
        message(paste("Was batch significant?",phenTestResult$model.effect.batch))
        message(paste("Was variance equal?",phenTestResult$model.effect.variance))
        
        if (phenTestResult$model.effect.interaction)
        {
            sexualDimorphism <- "yes"
        }
        else {
            sexualDimorphism <- "no"
        }
        message(paste("\nWas there evidence of sexual dimorphism? ",
                        sexualDimorphism," (p-value ",
                        sprintf("%.4f",round(phenTestResult$model.output.interaction,digits=4)),")",sep=""))
        
        if (!is.null(phenTestResult$model.output.percentageChanges)){
            if (phenTestResult$numberSexes==2){
                message(paste("Genotype percentage change Female:",
                                round(phenTestResult$model.output.percentageChanges[1],digits=2)),"%")
                message(paste("Genotype percentage change Male:",
                                round(phenTestResult$model.output.percentageChanges[2],digits=2)),"%")
            } 
            else {
              if ("Female" %in% levels(phenTestResult$model.dataset$Sex)){
                message(paste("Genotype percentage change Female:",
                              round(phenTestResult$model.output.percentageChanges[1],digits=2)),"%")
              }  
              else {
                message(paste("Genotype percentage change Male:",
                              round(phenTestResult$model.output.percentageChanges[1],digits=2)),"%")
                
              }
            }   
        }
        
        message(paste("\n",line,sep=""))
        message("Classification Tag")
        message(line)
        message(classificationTag(phenTestResult,phenotypeThreshold=phenotypeThreshold))
        
        message(paste("\n",line,sep=""))
        message("Model Output Summary")
        message(line)
        summary(phenTestResult$model.output)$tTable
    }
    
    else if (phenTestResult$method %in% c("FE")){
        colnum <- 1
        message(line)
        message(paste("Model Output ('*' highlights results with p-values less than threshold ",phenotypeThreshold
                        ,")",sep=""))
        message(line)
        #message("'All' data:")
        results_all <- matrix(0,2,1)
        results_all[1] <- sprintf("%.4f",round(phenTestResult$model.output$all$p.val,digits=4))
        #results_all[1] <- format(phenTestResult$model.output$all$p.val,  digits=4, scientific=F)
        results_all[2] <-  paste(format(phenTestResult$model.output$ES, nsmall = 0),"%",sep="")
        rownames(results_all) <- c("p-value:","Effect size:")
        colnames(results_all) <- c("")
        
        if (results_all[1] <= phenotypeThreshold){
            rownames(results_all)[1] <- paste("*",rownames(results_all)[1])
            rownames(results_all)[2] <- paste("*",rownames(results_all)[2])
        }
        else{
            rownames(results_all)[1] <- paste(" ",rownames(results_all)[1])
            rownames(results_all)[2] <- paste(" ",rownames(results_all)[2])
        }

        results <- results_all
        colnames(results)[colnum] <- "All"
        colnum <- colnum + 1
        #print(results_all)
        #message(paste("All p-val: ",phenTestResult$model.output$all$p.val,sep=""))
        
        #message(paste("All effect size: ",phenTestResult$model.output$ES,"%",sep=""))
        
        if (!is.null(phenTestResult$model.output$male)){
            #message("\n'Males only' data:")
            results_male <- matrix(0,2,1)
            results_male[1] <- sprintf("%.4f",round(phenTestResult$model.output$male$p.val, digits=4))
            results_male[2] <-  paste(format(phenTestResult$model.output$ES_male, nsmfemale = 0),"%",sep="")
            rownames(results_male) <- c("p-value:","Effect size:")
            colnames(results_male) <- c("")
            
            if (results_male[1] <= phenotypeThreshold){
                rownames(results_male)[1] <- paste("*",rownames(results_male)[1])
                rownames(results_male)[2] <- paste("*",rownames(results_male)[2])
            }
            else{
                rownames(results_male)[1] <- paste(" ",rownames(results_male)[1])
                rownames(results_male)[2] <- paste(" ",rownames(results_male)[2])
            }
            
            results <- cbind(results,results_male)
            colnames(results)[colnum] <- "Males only"
            colnum <- colnum + 1
            #print(results_male)
            #message(paste("Males only p-val: ",phenTestResult$model.output$male$p.val,sep=""))
            
            #message(paste("Males only effect size: ",phenTestResult$model.output$ES_male,"%",sep=""))
       
        }
        if (!is.null(phenTestResult$model.output$female)){
            #message("\n'Females only' data:")
            results_female <- matrix(0,2,1)
            results_female[1] <- sprintf("%.4f",round(phenTestResult$model.output$male$p.val, digits=4))
            results_female[2] <-  paste(format(phenTestResult$model.output$ES_male, nsmfemale = 0),"%",sep="")
            rownames(results_female) <- c("p-value:","Effect size:")
            colnames(results_female) <- c("")
            
            if (results_female[1] <= phenotypeThreshold){
                rownames(results_female)[1] <- paste("*",rownames(results_female)[1])
                rownames(results_female)[2] <- paste("*",rownames(results_female)[2])
            }
            else{
                rownames(results_female)[1] <- paste(" ",rownames(results_female)[1])
                rownames(results_female)[2] <- paste(" ",rownames(results_female)[2])
            }
            
            results <- cbind(results,results_female)
            colnames(results)[colnum] <- "Females only"
            colnum <- colnum + 1
            #print(results_female)
            #message(paste("Females only p-val: ",phenTestResult$model.output$female$p.val,sep=""))
            
            #message(paste("Females only effect size: ",phenTestResult$model.output$ES_female,"%",sep=""))                     
        }
        print(results,quote=FALSE)
        message(paste("\n",line,sep=""))
        message("Classification Tag")
        message(line)
        
        message(classificationTag(phenTestResult,phenotypeThreshold=phenotypeThreshold))

        message(paste("\n",line,sep=""))
        message("Count Matrices")
        message(line)
        message("\n'All' matrix:")
        print(phenTestResult$model.output$count_matrix_all)

        #message("\n'All' percentage matrix:")
        #print(phenTestResult$model.output$percentage_matrix_all)
        #message("\nMatrix 'all' statistics:")
        #print(phenTestResult$model.output$stat_all)
        
        if (!is.null(phenTestResult$model.output$male)){
            message("\n'Males only' matrix:")
            print(phenTestResult$model.output$count_matrix_male)
            #message("\n'Males only' percentage matrix:")
            #print(phenTestResult$model.output$percentage_matrix_male)
            #message("\nMatrix 'males only' statistics:")
            #print(phenTestResult$model.output$stat_male)  
            
        }
        if (!is.null(phenTestResult$model.output$female)){
            message("\n'Females only' matrix:")
            print(phenTestResult$model.output$count_matrix_female)
            #message("\n'Females only' percentage matrix:")
            #print(phenTestResult$model.output$percentage_matrix_female)
            #message("\nMatrix 'females only' statistics:")
            #print(phenTestResult$model.output$stat_female) 
        }
    }    
    
    else if (phenTestResult$method %in% c("RR")){
        colnum <- 1
        message(line)
        message(paste("Model Output ('*' highlights results with p-values less than threshold ",phenotypeThreshold
        ,")",sep=""))
        message(line)
        #message("'All' data:")
        
        RROutput <- phenTestResult$model.output$all
        if (RROutput[1] <= phenotypeThreshold){
            rownames(RROutput)[1] <- paste("*",rownames(RROutput)[1])
            rownames(RROutput)[2] <- paste("*",rownames(RROutput)[2])
            #RROutput[1] <- paste("*",RROutput[1])
            #RROutput[2] <- paste("*",RROutput[2])
        }
        else {
            rownames(RROutput)[1] <- paste(" ",rownames(RROutput)[1])
            rownames(RROutput)[2] <- paste(" ",rownames(RROutput)[2])
        }
        if (RROutput[3] <= phenotypeThreshold){
            rownames(RROutput)[3] <- paste("*",rownames(RROutput)[3])
            rownames(RROutput)[4] <- paste("*",rownames(RROutput)[4])
            #RROutput[3] <- paste("*",RROutput[3])
            #RROutput[4] <- paste("*",RROutput[4])
        }
        else{
            rownames(RROutput)[3] <- paste(" ",rownames(RROutput)[3])
            rownames(RROutput)[4] <- paste(" ",rownames(RROutput)[4])
        }
        #print(RROutput)

        RROutput2 <- RROutput
        colnames(RROutput2)[colnum] <- "All"
        colnum <- colnum + 1
        
        if (!is.null(phenTestResult$model.output$male)){
            #message("\n'Males only' data:")
            RROutput <- phenTestResult$model.output$male
            if (RROutput[1] <= phenotypeThreshold){
                rownames(RROutput)[1] <- paste("*",rownames(RROutput)[1])
                rownames(RROutput)[2] <- paste("*",rownames(RROutput)[2])
                #RROutput[1] <- paste("*",RROutput[1])
                #RROutput[2] <- paste("*",RROutput[2])
            }
            else{
                rownames(RROutput)[1] <- paste(" ",rownames(RROutput)[1])
                rownames(RROutput)[2] <- paste(" ",rownames(RROutput)[2])
            }
            if (RROutput[3] <= phenotypeThreshold){
                rownames(RROutput)[3] <- paste("*",rownames(RROutput)[3])
                rownames(RROutput)[4] <- paste("*",rownames(RROutput)[4])
                #RROutput[3] <- paste("*",RROutput[3])
                #RROutput[4] <- paste("*",RROutput[4])
            }
            else {
                rownames(RROutput)[3] <- paste(" ",rownames(RROutput)[3])
                rownames(RROutput)[4] <- paste(" ",rownames(RROutput)[4])
            }
            #print(RROutput)
            RROutput2 <- cbind(RROutput2,RROutput)
            colnames(RROutput2)[colnum] <- "Males only"
            colnum <- colnum + 1
        }
        if (!is.null(phenTestResult$model.output$female)){
            #message("\n'Females only' data:")
            RROutput <- phenTestResult$model.output$female
            if (RROutput[1] <= phenotypeThreshold){
                #rownames(RROutput)[1] <- paste("*",rownames(RROutput)[1])
                #rownames(RROutput)[2] <- paste("*",rownames(RROutput)[2])
                rownames(RROutput)[1] <- paste("*",rownames(RROutput)[1])
                rownames(RROutput)[2] <- paste("*",rownames(RROutput)[2])
            }
            else{
                rownames(RROutput)[1] <- paste(" ",rownames(RROutput)[1])
                rownames(RROutput)[2] <- paste(" ",rownames(RROutput)[2])
            }
            if (RROutput[3] <= phenotypeThreshold){
                #rownames(RROutput)[3] <- paste("*",rownames(RROutput)[3])
                #rownames(RROutput)[4] <- paste("*",rownames(RROutput)[4])
                rownames(RROutput)[3] <- paste("*",rownames(RROutput)[3])
                rownames(RROutput)[4] <- paste("*",rownames(RROutput)[4])
            }
            else{
                rownames(RROutput)[3] <- paste(" ",rownames(RROutput)[3])
                rownames(RROutput)[4] <- paste(" ",rownames(RROutput)[4])
            }
            #print(RROutput)
            RROutput2 <- cbind(RROutput2,RROutput)
            colnames(RROutput2)[colnum] <- "Females only"
            colnum <- colnum + 1
        }

        print(RROutput2,quote=FALSE)
        message(paste("\n",line,sep=""))
        message("Classification Tag")
        message(line)
        message(classificationTag(phenTestResult,phenotypeThreshold=phenotypeThreshold))
        
        message(paste("\n",line,sep=""))
        message("Thresholds")
        message(line)
        thresholds <- phenTestResult$model.output.quality
        thresholds <- rbind(phenotypeThreshold,thresholds)
        rownames(thresholds)[1] <- "p-value threshold:"
        print(thresholds,quote=FALSE)
        message(paste("\n",line,sep=""))
        message("Count Matrices")
        message(line)
        message("'All' matrix:")
        print(phenTestResult$model.output$count_matrix_all)        
        if (!is.null(phenTestResult$model.output$male)){
            message("\n'Males only' matrix:")
            print(phenTestResult$model.output$count_matrix_male)
        }
        if (!is.null(phenTestResult$model.output$female)){
            message("\n'Females only' matrix:")
            print(phenTestResult$model.output$count_matrix_female)
        }

    }  
}

##------------------------------------------------------------------------------
## Generate graphs and store them in the defined directory
generateGraphs <- function(phenTestResult, dir, graphingName=NULL, type="Xlib")
{
    if (phenTestResult$method=="MM"){
        
        if (is.null(graphingName))
        graphingName = phenTestResult$depVariable
        
        
        #1
        graph_name=file.path(dir, "boxplotSexGenotype.png")
        png(graph_name,type=type)
        boxplotSexGenotype(phenTestResult$model.dataset,
                phenTestResult$depVariable, graphingName=graphingName)
        
        dev.off()
        
        #2
        if (('Batch' %in% colnames(phenTestResult$model.dataset))){
            graph_name=file.path(dir, "boxplotSexGenotypeBatch.png")
            png(graph_name,type=type)
            boxplotSexGenotypeBatch(phenTestResult$model.dataset,
                    phenTestResult$depVariable, graphingName=graphingName)
            
            dev.off()
        }
        
        #3
        if (('Weight' %in% colnames(phenTestResult$model.dataset))){
            graph_name=file.path(dir, "scatterplotGenotypeWeight.png")
            png(graph_name,type=type)
            scatterplotGenotypeWeight(phenTestResult$model.dataset,
                    phenTestResult$depVariable, graphingName=graphingName)
            
            dev.off()
        }
        
        
        
        #4
        graph_name=file.path(dir, "qqplotGenotype.png")
        png(graph_name,type=type)
        qqplotGenotype(phenTestResult,outputMessages=FALSE)
        dev.off()
        
        #5
        graph_name=file.path(dir, "plotResidualPredicted.png")
        png(graph_name,type=type)
        plotResidualPredicted(phenTestResult)
        dev.off()
        
        if (('Batch' %in% colnames(phenTestResult$model.dataset)) 
                && phenTestResult$model.effect.batch){
            graph_name=file.path(dir, "qqplotRandomEffects.png")
            png(graph_name,type=type)
            qqplotRandomEffects(phenTestResult,outputMessages=FALSE)
            
            dev.off()
        }    
        
        #7
        if (('Batch' %in% colnames(phenTestResult$model.dataset))){
            graph_name=file.path(dir, "boxplotResidualBatch.png")
            png(graph_name,type=type)
            boxplotResidualBatch(phenTestResult,outputMessages=FALSE)
            
            dev.off()
        }    
        
        #8
        if (('Batch' %in% colnames(phenTestResult$model.dataset)) 
                && phenTestResult$model.effect.batch){
            graph_name=file.path(dir, "qqplotRotatedResiduals.png")
            png(graph_name,type=type)
            qqplotRotatedResiduals(phenTestResult,outputMessages=FALSE)
            
            dev.off()
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