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
# simpleOutput.R contains simpleOutput, generateGraphs functions
#-----------------------------------------------------------------------------------
summaryOutput <- function(phenTestResult,phenotypeThreshold=0.01)
# Wrapper to prepare the output of the modelling and testing results in simple user friendly form. Assumes that modelling results are 
# stored in the phenTestResult object (output from functions testDataset and buildFinalModel)
{
    message(paste("Test for dependent variable: ",phenTestResult$depVariable,sep=""))
    
    message(paste("Method: ",switch(phenTestResult$method,MM = "Mixed Model framework",FE = "Fisher Exact Test framework"),"\n",sep=""))

    if (phenTestResult$method=="MM") {
        message(paste("Was batch significant?",phenTestResult$model.effect.batch))
        
        message(paste("Was variance equal?",phenTestResult$model.effect.variance))
        
        if (phenTestResult$model.effect.interaction)
            sexualDimorphism = "yes"
        else 
            sexualDimorphism = "no"
        message(paste("Was there evidence of sexual dimorphism? ",sexualDimorphism," (p-value ",round(phenTestResult$model.output.interaction,digits=3),")",sep=""))
        
        message(paste("Final fitted model:",format(phenTestResult$model.formula.genotype)))
        
        message("Model output:")
        
        message(paste("Genotype effect:",round(phenTestResult$model.output.genotype.nulltest.pVal,digits=9)))
        
        message(paste("Classification tag:", classificationTag(phenTestResult,phenotypeThreshold=phenotypeThreshold)))
        
        summary(phenTestResult$model.output)$tTable
    }
    
    else {
        message("Model output:")
        
        message(paste("All data p-val: ",round(phenTestResult$model.output$all$p.val,digits=3),sep=""))
        
        message(paste("All data effect size: ",round(phenTestResult$model.output$stat_all$cramer,digits=3),sep=""))
        
        if (!is.null(phenTestResult$model.output$male)){
            message(paste("Males only p-val: ",round(phenTestResult$model.output$male$p.val,digits=3),sep=""))
            message(paste("Males only effect size: ",round(phenTestResult$model.output$stat_male$cramer,digits=3),sep=""))
            
        }
        if (!is.null(phenTestResult$model.output$female)){
            message(paste("Females only p-val: ",round(phenTestResult$model.output$female$p.val,digits=3),sep=""))
            message(paste("Females only effect size: ",round(phenTestResult$model.output$stat_female$cramer,digits=3),sep=""))
           
        }
        
    }
    
    
    
    
}
#-----------------------------------------------------------------------------------
generateGraphs <- function(phenList, phenTestResult, dir, graphingName=NULL, type="Xlib")
# Generate all possible graphs and store them in the defined directory
{
    # ADD graphing name and call all graphs 
    if (is.null(graphingName))
        graphingName = phenTestResult$depVariable

    
    #1
    graph_name=file.path(dir, "boxplotGenderGenotype.png")
    png(graph_name,type=type)
    boxplotGenderGenotype(phenList,phenTestResult$depVariable, graphingName=graphingName)
    dev.off()
    
    #2
    if (('Batch' %in% colnames(phenList$dataset))){
        graph_name=file.path(dir, "boxplotGenderGenotypeBatch.png")
        png(graph_name,type=type)
        boxplotGenderGenotypeBatch(phenList,phenTestResult$depVariable, graphingName=graphingName)
        dev.off()
    }
    
    #3
    if (('Weight' %in% colnames(phenList$dataset))){
        graph_name=file.path(dir, "scatterplotGenotypeWeight.png")
        png(graph_name,type=type)
        scatterplotGenotypeWeight(phenList,phenTestResult$depVariable, graphingName=graphingName)
        dev.off()
    }
    
    if (phenTestResult$method=="MM"){
    
        #4
        graph_name=file.path(dir, "qqplotGenotype.png")
        png(graph_name,type=type)
        qqplotGenotype(phenTestResult)
        dev.off()
        
        #5
        graph_name=file.path(dir, "plotResidualPredicted.png")
        png(graph_name,type=type)
        plotResidualPredicted(phenTestResult)
        dev.off()
        
        if (('Batch' %in% colnames(phenList$dataset)) && phenTestResult$model.effect.batch){
            graph_name=file.path(dir, "qqplotRandomEffects.png")
            png(graph_name,type=type)
            qqplotRandomEffects(phenTestResult)
            dev.off()
        }    
        
        #7
        if (('Batch' %in% colnames(phenList$dataset))){
            graph_name=file.path(dir, "boxplotResidualBatch.png")
            png(graph_name,type=type)
            boxplotResidualBatch(phenTestResult)
            dev.off()
        }    
        
        #8
        if (('Batch' %in% colnames(phenList$dataset)) && phenTestResult$model.effect.batch){
            graph_name=file.path(dir, "qqplotRotatedResiduals.png")
            png(graph_name,type=type)
            qqplotRotatedResiduals(phenTestResult)
            dev.off()
        }   
        
    } 
    
    
}    
 