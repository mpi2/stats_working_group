# graphsDataset.R contains functions for dataset diagnostic plots:
# boxplotGenderGenotype, boxplotGenderGenotypeBatch, scatterplotGenotypeWeigh

#-----------------------------------------------------------------------------------
# Raw data boxplot split by gender and genotype
boxplotGenderGenotype<-function(phenList, depVariable){
    if(is(phenList,"PhenList")) {
        x <- phenList$phendata       
        
    } else {
        x <- as.data.frame(phenList)
    }
    if (is.null(depVariable)) stop("Please define dependant variable 'depVariable'")
    if (!is(x$Genotype)) stop("Genotype values are not defined")  
    if (!is(x$Gender)) stop("Gender values are not defined")    
    
    numberofgenders=length(levels(x$Gender))
    if(numberofgenders==2){
        Male = subset(x, x$Gender=="Male")
        Female= subset(x, x$Gender=="Female")
        par(mfrow = c(1, 2))
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
# Row data boxplot by gender,genotype and batch (Assay.date)
boxplotGenderGenotypeBatch<-function(phenList, depVariable){
    if(is(phenList,"PhenList")) {
        x <- phenList$phendata
        refGenotype <- phenList$refGenotype
        
    } else {
        x <- as.data.frame(phenList)
        refGenotype <- "+/+"
    }
    if (is.null(depVariable)) stop("Please define dependant variable 'depVariable'")
    if (!is(x$Genotype)) stop("Genotype values are not defined")  
    if (!is(x$Gender)) stop("Gender values are not defined") 
    
    numberofgenders=length(levels(x$Gender))
    y_range=c(min(x[ ,depVariable], na.rm=TRUE), max((x[ ,depVariable]), na.rm=TRUE))
    x_range=c(levels(x$Assay.Date))

    if(numberofgenders==2){
        Male = subset(x, x$Gender=="Male")
        Female= subset(x, x$Gender=="Female")
        par(mfrow = c(1, 2))
        boxplot(Male[ , depVariable]~Male$Genotype+Male$Assay.Date, subset=(Male$Genotype=="+/+"), ylab=depVariable, ylim=y_range, xlab="Genotype", names=NULL)
        boxplot(Male[ , depVariable]~Male$Genotype + Male$Assay.Date, add=TRUE, subset=(Male$Genotype!="+/+"), ylim=y_range, ylab=depVariable, xlab="Genotype", col="red", names=NULL)
        legend("topright", "Male", cex=1.3, bty="n")
        boxplot(Female[ , depVariable]~Female$Genotype + Female$Assay.Date, subset=(Female$Genotype=="+/+"),ylim=y_range,ylab=depVariable, xlab="Genotype", names=NULL )
        boxplot(Female[ , depVariable]~Female$Genotype + Female$Assay.Date, add=TRUE, subset=(Female$Genotype!="+/+"),ylim=y_range, ylab=depVariable, xlab="Genotype", col="red",names=NULL)
        legend("topright", "Female", cex=1.3, bty="n")
    }else{
        par(mfrow=c(1,1))
        boxplot(x[ ,depVariable]~x$Genotype+Male$Assay.Date,subset=(Male$Genotype=="+/+"), ylab=depVariable, xlab="Genotype", names=NULL) # xlim=x_range,
        boxplot(x[ ,depVariable]~x$Genotype+Male$Assay.Date,subset=(Male$Genotype!="+/+"), ylab=depVariable, xlab="Genotype", col="red", names=NULL)    
    }

}    

#-----------------------------------------------------------------------------------
# Body weight versus dependent variable scatter plot
scatterplotGenotypeWeight<-function(phenList, depVariable){
    require(car)
    if(is(phenList,"PhenList")) {
        x <- phenList$phendata     
        
    } else {
        x <- as.data.frame(phenList)
    }
    if (is.null(depVariable)) stop("Please define dependant variable 'depVariable'")
    
    model.formula <- as.formula(paste(depVariable, "~", paste("Weight", "Genotype", sep= "|")))
    scatterplot(data=x, model.formula)
}

#-----------------------------------------------------------------------------------




