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
## simpleOutput.R contains simpleOutput, generateGraphs functions
##------------------------------------------------------------------------------
## Wrapper to prepare the output of the modelling and testing results in simple 
## user friendly form. Assumes that modelling results are stored in the
## phenTestResult object (output from functions testDataset and buildFinalModel)
summaryOutput <- function(phenTestResult,phenotypeThreshold=0.01)
{
    message(paste("Test for dependent variable: ",
                    phenTestResult$depVariable,sep=""))
    
    message(paste("Method: ",switch(phenTestResult$method,
                            MM = "Mixed Model framework",FE = "Fisher Exact Test framework",
                            RR = "Reference Ranges Plus framework"),"\n",sep=""))
    
    if (phenTestResult$method=="MM") {
        message(paste("Was batch significant?",phenTestResult$model.effect.batch))
        
        message(paste("Was variance equal?",phenTestResult$model.effect.variance))
        
        if (phenTestResult$model.effect.interaction)
        sexualDimorphism = "yes"
        else 
        sexualDimorphism = "no"
        message(paste("Was there evidence of sexual dimorphism? ",
                        sexualDimorphism," (p-value ",
                        round(phenTestResult$model.output.interaction,digits=3),")",sep=""))
        
        message(paste("Final fitted model:",
                        format(phenTestResult$model.formula.genotype)))
        
        message("Model output:")
        
        message(paste("Genotype effect:",
                        round(phenTestResult$model.output.genotype.nulltest.pVal,digits=9)))
        
        message(paste("Classification tag:", 
                        classificationTag(phenTestResult,phenotypeThreshold=phenotypeThreshold)))
        
        summary(phenTestResult$model.output)$tTable
    }
    
    else if (phenTestResult$method %in% c("FE")){
        message("Model output:")
        
        message(paste("All data p-val: ",
                        phenTestResult$model.output$all$p.val,sep=""))
        
        message(paste("All data effect size: ",
                        phenTestResult$model.output$ES,"%",sep=""))
        
        if (!is.null(phenTestResult$model.output$male)){
            message(paste("Males only p-val: ",
                            phenTestResult$model.output$male$p.val,sep=""))
            
            message(paste("Males only effect size: ",
                            phenTestResult$model.output$ES_male,"%",sep=""))
        }
        if (!is.null(phenTestResult$model.output$female)){
            message(paste("Females only p-val: ",
                            phenTestResult$model.output$female$p.va,sep=""))
            
            message(paste("Females only effect size: ",
                            phenTestResult$model.output$ES_female,"%",sep=""))                     
        }
        message(paste("Classification tag:", 
                        classificationTag(phenTestResult)))
        
        ## Matrices and statistics
        message("\nMatrix 'all':")
        print(phenTestResult$model.output$count_matrix_all)
        message("\nPercentage matrix 'all' statistics:")
        print(phenTestResult$model.output$percentage_matrix_all)
        message("\nMatrix 'all' statistics:")
        print(phenTestResult$model.output$stat_all)
        
        if (!is.null(phenTestResult$model.output$male)){
            message("\nMatrix 'males only':")
            print(phenTestResult$model.output$count_matrix_male)
            message("\nPercentage matrix 'males only' statistics:")
            print(phenTestResult$model.output$percentage_matrix_male)
            message("\nMatrix 'males only' statistics:")
            print(phenTestResult$model.output$stat_male)  
            
        }
        if (!is.null(phenTestResult$model.output$female)){
            message("\nMatrix 'females only':")
            print(phenTestResult$model.output$count_matrix_female)
            message("\nPercentage matrix 'females only' statistics:")
            print(phenTestResult$model.output$percentage_matrix_female)
            message("\nMatrix 'females only' statistics:")
            print(phenTestResult$model.output$stat_female) 
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