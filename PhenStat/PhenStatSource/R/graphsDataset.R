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
## graphsDataset.R contains functions for dataset diagnostic plots:
## boxplotGenderGenotype, boxplotGenderGenotypeBatch, 
## scatterplotGenotypeWeigh & categoricalBarplot
##------------------------------------------------------------------------------
## Raw data boxplot: split by gender and genotype
boxplotGenderGenotype<-function(phenList, depVariable=NULL, graphingName=NULL){
    
    ## Checks
    if(is(phenList,"PhenList")) {
        x <- phenList$dataset               
    } else {
        x <- phenList
    }
    if (is.null(depVariable)) 
        stop("Please define dependent variable 'depVariable'.")
    
    if (is.null(graphingName))
        graphingName <- depVariable
    
    if (!(depVariable %in% colnames(x)))
    stop(paste(depVariable,"column is missed in the dataset."))
    else {
        columnOfInterest <- x[,c(depVariable)]
        
        ## Test: depVariable is numeric 
        if(!is.numeric(columnOfInterest))
        stop(paste(depVariable,"variable is not numeric. 
                        Can't create a plot based on it."))
    }       
    
    
    ## Plot creation
    numberofgenders <- length(levels(x$Gender))
    if(numberofgenders==2){
        Male <- subset(x, x$Gender=="Male")
        Female <- subset(x, x$Gender=="Female")
        op <- par(mfrow=c(1,2))
        boxplot(Male[ , depVariable]~Male$Genotype, 
                ylab=graphingName, xlab="Genotype")
        legend("topright", "Male", cex=1.3, bty="n")
        boxplot(Female[ , depVariable]~Female$Genotype, 
                ylab=graphingName, xlab="Genotype")
        legend("topright", "Female", cex=1.3, bty="n")
        par(op) 
        op_normal <- par(mfrow=c(1,1))
        par(op_normal) 
    }else{
        op <- par(mfrow=c(1,1))
        boxplot(x[ ,depVariable]~x$Genotype, ylab=graphingName, xlab="Genotype")  
        par(op)  
    }
}    
##------------------------------------------------------------------------------
## Raw data boxplot: split by gender,genotype and batch 
boxplotGenderGenotypeBatch<-function(phenList, depVariable=NULL, graphingName=NULL){
    
    ## Checks
    if(is(phenList,"PhenList")) {
        x <- phenList$dataset
        refGenotype <- phenList$refGenotype        
    } else {
        x <- phenList
    }

            
    if (is.null(depVariable)) 
            stop("Please define dependent variable 'depVariable'.")
    
    if (is.null(graphingName))
            graphingName <- depVariable
    
    
    if (!(depVariable %in% colnames(x)))
    stop(paste(depVariable,"column is missed in the dataset."))
    else {
        columnOfInterest <- x[,c(depVariable)]
        
        ## Test: depVariable is numeric 
        if(!is.numeric(columnOfInterest))
        stop(paste(depVariable,"variable is not numeric. 
                        Can't create a plot based on it."))
    }       
    
    
    if (!('Batch' %in% colnames(x)))
    stop(paste("Batch column is missed in the dataset."))
    
    ## Plot creation
    numberofgenders <- length(levels(x$Gender))
    if (is.numeric(x[ ,depVariable]))   
    y_range <- c(min(x[ ,depVariable], na.rm=TRUE), 
            max((x[ ,depVariable]), na.rm=TRUE))
    else
    y_range <- c(1, length(levels(x[ ,depVariable])))
    
    if(numberofgenders==2){
        Male <- subset(x, x$Gender=="Male")
        Female <- subset(x, x$Gender=="Female")
        Male$Batch <- factor(Male$Batch)
        Female$Batch <- factor(Female$Batch)
        
        op <- par(mfrow=c(1,2)) 
        
        boxplot(Male[ , depVariable]~Male$Genotype+Male$Batch, 
                subset=(Male$Genotype==phenList$refGenotype), 
                ylab=graphingName, ylim=y_range, xlab="Batch",  names=NULL)
        
        boxplot(Male[ , depVariable]~Male$Genotype + Male$Batch, add=TRUE, 
                subset=(Male$Genotype!=phenList$refGenotype), ylim=y_range, 
                ylab=graphingName, xlab="Batch",  col="red", names=NULL)
        
        legend("topright", "Male", cex=1.3, bty="n")
        
        boxplot(Female[ , depVariable]~Female$Genotype + Female$Batch, 
                subset=(Female$Genotype==phenList$refGenotype),ylim=y_range,
                ylab=graphingName, xlab="Batch",  names=NULL )
        
        boxplot(Female[ , depVariable]~Female$Genotype + Female$Batch, add=TRUE, 
                subset=(Female$Genotype!=phenList$refGenotype),ylim=y_range, 
                ylab=graphingName, xlab="Batch",  col="red",names=NULL)
        
        legend("topright", "Female", cex=1.3, bty="n")
        
        par(op)
        
        op_normal <- par(mfrow=c(1,1))
        par(op_normal) 
    }else{
        op <- par(mfrow=c(1,1))
        boxplot(x[ ,depVariable]~x$Genotype+x$Batch,
                subset=(x$Genotype==phenList$refGenotype), 
                ylab=graphingName, xlab="Batch", names=NULL) 
        
        boxplot(x[ ,depVariable]~x$Genotype+x$Batch,
                subset=(x$Genotype!=phenList$refGenotype), add=TRUE, 
                ylab=graphingName, xlab="Batch",  col="red", names=NULL)    
        
        par(op)
    }
    
}    

##------------------------------------------------------------------------------
## Raw data scatterplot: body weight versus dependant variable
scatterplotGenotypeWeight<-function(phenList, depVariable=NULL, graphingName=NULL){
    
    ## Checks
    if(is(phenList,"PhenList")) {
        x <- phenList$dataset     
        
    } else {
        x <- phenList
    }
    if (is.null(depVariable)) 
        stop("Please define dependent variable 'depVariable'.")
    
    if (is.null(graphingName))
        graphingName <- depVariable
    
    
    if (!(depVariable %in% colnames(x)))
    stop(paste(depVariable,"column is missed in the dataset."))
    else {
        columnOfInterest <- x[,c(depVariable)]        
        ## Test: depVariable is numeric 
        if(!is.numeric(columnOfInterest))
        stop(paste(depVariable,"variable is not numeric. 
                        Can't create a plot based on it."))
    }       
    
    
    if (!('Weight' %in% colnames(x)))
    stop("Weight is missed in the dataset.")
    
    ## Plot creation
    model.formula <- as.formula(paste(depVariable, "~", 
                    paste("Weight", "Genotype", sep= "|")))
    
    scatterplot(data=x, model.formula, ylab=graphingName)
}
##------------------------------------------------------------------------------