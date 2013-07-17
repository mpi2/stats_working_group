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

# graphsDataset.R contains functions for dataset diagnostic plots:
# boxplotGenderGenotype, boxplotGenderGenotypeBatch, scatterplotGenotypeWeigh

#-----------------------------------------------------------------------------------
# Raw data boxplot: split by gender and genotype
boxplotGenderGenotype<-function(phenList, depVariable){
    if(is(phenList,"PhenList")) {
        x <- phenList$dataset       
        
    } else {
        stop("Please create a PhenList object first")
    }
    if (is.null(depVariable)) stop("Please define dependent variable 'depVariable'")
    
    numberofgenders=length(levels(x$Gender))
    if(numberofgenders==2){
        Male = subset(x, x$Gender=="Male")
        Female= subset(x, x$Gender=="Female")
        par(mfrow=c(1,2)) 
        boxplot(Male[ , depVariable]~Male$Genotype, ylab=depVariable, xlab="Genotype")
        legend("topright", "Male", cex=1.3, bty="n")
        boxplot(Female[ , depVariable]~Female$Genotype, ylab=depVariable, xlab="Genotype")
        legend("topright", "Female", cex=1.3, bty="n")
    }else{
        par(mfrow=c(1,1))
        boxplot(x[ ,depVariable]~x$Genotype, ylab=depVariable, xlab="Genotype")    
    }
}    
#-----------------------------------------------------------------------------------
# Raw data boxplot: split by gender,genotype and batch 
boxplotGenderGenotypeBatch<-function(phenList, depVariable){
    if(is(phenList,"PhenList")) {
        x <- phenList$dataset
        refGenotype <- phenList$refGenotype
        
    } else {
        stop("Please create a PhenList object first")
    }
    if (is.null(depVariable)) stop("Please define dependent variable 'depVariable'")
    
    numberofgenders=length(levels(x$Gender))
    y_range=c(min(x[ ,depVariable], na.rm=TRUE), max((x[ ,depVariable]), na.rm=TRUE))
    x_range=c(levels(x$Batch))
    
    if(numberofgenders==2){
        Male = subset(x, x$Gender=="Male")
        Female= subset(x, x$Gender=="Female")
        par(mfrow=c(1,2)) 
        boxplot(Male[ , depVariable]~Male$Genotype+Male$Batch, subset=(Male$Genotype=="+/+"), ylab=depVariable, ylim=y_range, xlab="Genotype", names=NULL)
        boxplot(Male[ , depVariable]~Male$Genotype + Male$Batch, add=TRUE, subset=(Male$Genotype!="+/+"), ylim=y_range, ylab=depVariable, xlab="Genotype", col="red", names=NULL)
        legend("topright", "Male", cex=1.3, bty="n")
        boxplot(Female[ , depVariable]~Female$Genotype + Female$Batch, subset=(Female$Genotype=="+/+"),ylim=y_range,ylab=depVariable, xlab="Genotype", names=NULL )
        boxplot(Female[ , depVariable]~Female$Genotype + Female$Batch, add=TRUE, subset=(Female$Genotype!="+/+"),ylim=y_range, ylab=depVariable, xlab="Genotype", col="red",names=NULL)
        legend("topright", "Female", cex=1.3, bty="n")
    }else{
        par(mfrow=c(1,1))
        boxplot(x[ ,depVariable]~x$Genotype+x$Batch,subset=(x$Genotype=="+/+"), ylab=depVariable, xlab="Genotype", names=NULL) # xlim=x_range,
        boxplot(x[ ,depVariable]~x$Genotype+x$Batch,subset=(x$Genotype!="+/+"), ylab=depVariable, xlab="Genotype", col="red", names=NULL)    
    }
    
}    

#-----------------------------------------------------------------------------------
# Raw data scatterplot: body weight versus dependant variable
scatterplotGenotypeWeight<-function(phenList, depVariable){
    require(car)
    if(is(phenList,"PhenList")) {
        x <- phenList$dataset     
        
    } else {
        stop("Please create a PhenList object first")
    }
    if (is.null(depVariable)) stop("Please define dependent variable 'depVariable'")
    
    model.formula <- as.formula(paste(depVariable, "~", paste("Weight", "Genotype", sep= "|")))
    scatterplot(data=x, model.formula)
}


#-----------------------------------------------------------------------------------




