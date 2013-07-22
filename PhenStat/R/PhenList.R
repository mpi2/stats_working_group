# Copyright ¬© 2011-2013 EMBL - European Bioinformatics Institute
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


# PhenList.R contains PhenList and checkDataset functions to construct new PhenList object from components and to check 
# dataset integrity

PhenList <- function(dataset=matrix(0,0,0), outputMessages=TRUE, refGenotype='+/+', dataset.stat=NULL, dataset.colname.batch=NULL,
dataset.colname.genotype=NULL, dataset.colname.gender=NULL, dataset.colname.weight=NULL, dataset.values.missingValue=NULL, 
dataset.values.male=NULL, dataset.values.female=NULL, dataset.clean=TRUE, testGenotype=NULL, hemiGenotype=NULL) 

# Construct PhenList object from components with data quality checks

{
    dataset <- dataset[,order(names(dataset))]
    
    
# Rename columns if needed  
    if (dataset.clean){  
        if(!is.null(dataset.colname.batch)) colnames(dataset)[colnames(dataset) == dataset.colname.batch] <-'Batch'
        else { 
            if (!('Assay.Date' %in% colnames(dataset))){
                colnames(dataset)[colnames(dataset) == 'Assay.Date'] <-'Batch'
                if (outputMessages)
                message("Warning: Dataset's column 'Assay.Date' has been renamed to 'Batch' and will be used for the Batch effect modeling.")
            }
            else          
			
			if (length(colnames(dataset)[grep("batch", tolower(colnames(dataset)))])>0 && outputMessages){
				batch_potential_columns<-paste(colnames(dataset)[grep("batch", tolower(colnames(dataset)))], collapse="', '" )  
				message(paste("Warning: Dataset contains columns that might be used for Batch effect modeling, for instance '",batch_potential_columns,"'.",sep=""))
			}
			
        }
        if(!is.null(dataset.colname.genotype)) colnames(dataset)[colnames(dataset) == dataset.colname.genotype] <-'Genotype'
        if(!is.null(dataset.colname.gender)) colnames(dataset)[colnames(dataset) == dataset.colname.gender] <-'Gender'
        if(!is.null(dataset.colname.weight)) colnames(dataset)[colnames(dataset) == dataset.colname.weight] <-'Weight'
        
# Replace missing values specified in the user format with NA if needed 
        if(!is.null(dataset.values.missingValue)) dataset[dataset == dataset.values.missingValue] <- NA 
    }
    
# Clean the empty records   
    dataset<-dataset[dataset$Gender!="",]
    dataset<-dataset[dataset$Genotype!="",]
    if ('Batch' %in% colnames(dataset))
	dataset<-dataset[dataset$Batch!="",]
    
# Renew levels 
    dataset$Gender<-factor(dataset$Gender)
    dataset$Genotype<-factor(dataset$Genotype)
    if ('Batch' %in% colnames(dataset))
	dataset$Batch<-factor(dataset$Batch)
	
    
# Replace values for genders with 'Male','Female' if needed 
    if (dataset.clean){
        if(!is.null(dataset.values.female)) levels(dataset$Gender)[levels(dataset$Gender)==dataset.values.female] <- "Female"
        if(!is.null(dataset.values.male)) levels(dataset$Gender)[levels(dataset$Gender)==dataset.values.male] <- "Male"
        
# Hemi to test replacement
        if (!is.null(hemiGenotype) && !is.null(testGenotype)) { 
            levels(dataset$Genotype)[levels(dataset$Genotype)==hemiGenotype] <- testGenotype
            if (outputMessages)
            message(paste("Warning: Hemizygotes '",hemiGenotype,"' have been relabled to homozygotes '",testGenotype,"'. If you don't want this behaviour please delete the hemizygotes records from the dataset.",sep=""))
        }
        
# Clean genotypes
        if (!is.null(testGenotype)){
            dataset <- dataset[-(dataset$Genotype!=testGenotype || dataset$Genotype!=refGenotype),]
            if (outputMessages)
			message(paste("Warning: Dataset has been cleaned, filtered out rows with genotype value other than ", testGenotype,"or",refGenotype))
        }
        
#Renew levels
        dataset$Genotype<-factor(dataset$Genotype)
        dataset$Gender<-factor(dataset$Gender)
    }    
	
	
    dataset <- checkDataset(dataset, outputMessages, refGenotype, dataset.clean)
    
    Genotype_levels=levels(dataset$Genotype)
    Gender_levels=levels(dataset$Gender)
    
# Statistics
    dataset.stat <- data.frame(Variables = colnames(dataset),Numeric = sapply(dataset, is.numeric), 
							   Continuous = sapply(dataset, function(x) if(is.numeric(x)) {if (length(unique(x))/length(x)>0.05) TRUE else FALSE} else FALSE),
							   Levels = sapply(dataset, function(x) if (length(unique(x))<4) paste(levels(x),collapse="*") else length(unique(x)) ), 
#UniqueObs = sapply(dataset, function(x) if(is.numeric(x)) {paste(round(length(unique(x))/length(x),digits=2)*100,"%",sep="")} else NA), 
							   NObs = sapply(dataset, function(x) length(na.omit(x))),
							   Mean = sapply(dataset, function(x) if(is.numeric(x)) round(mean(na.omit(x)),digits=2) else NA),
							   StdDev = sapply(dataset, function(x) if(is.numeric(x)) round(sd(na.omit(x)),digits=2) else NA),
							   Minimum = sapply(dataset, function(x) if(is.numeric(x)) round(min(na.omit(x)),digits=2) else NA),
							   Maximum = sapply(dataset, function(x) if(is.numeric(x)) round(max(na.omit(x)),digits=2) else NA))
    rownames(dataset.stat) <- NULL
    
	
	
    x <- new("PhenList",list(dataset=dataset))
    
    x$refGenotype <- refGenotype
    x$testGenotype <- testGenotype
    x$hemiGenotype <- hemiGenotype
    x$dataset.stat <- dataset.stat
    x$dataset.colname.batch <- dataset.colname.batch
    x$dataset.colname.genotype <- dataset.colname.genotype
    x$dataset.colname.gender <- dataset.colname.gender
    x$dataset.colname.weight <- dataset.colname.weight
    x$dataset.values.missingValue <- dataset.values.missingValue
    x$dataset.values.male <- dataset.values.male
    x$dataset.values.female <- dataset.values.female
    x$dataset.clean <- dataset.clean
    
    
    x
}

checkDataset <- function(dataset, outputMessages=TRUE, refGenotype="+/+", dataset.clean=FALSE)

# Check dataset for the minimum required info

{
    message <- ""
    message_dp <- ""
    pass <- TRUE
    
    
    nvar <- ncol(dataset)
    ntags <- nrow(dataset)
	
# Column names should be given
    if(nvar>0 && is.null(colnames(dataset))) {
        pass <- FALSE
        message <- paste(message, "Dataset with no column names.\n",sep="")
    }    
	
    
# Minimum required data
#if (!('Batch' %in% colnames(dataset))){
#    pass <- FALSE
#    message <- paste(message, "Dataset's 'Batch' column is missed\n")
#}    
    
# Check for mandatory columns: Genotype and Gender
    if (!('Genotype' %in% colnames(dataset))) {
        pass <- FALSE
        message <- paste(message,"Dataset's 'Genotype' column is missed.\n",sep="")
    }
	
    if (!('Gender' %in% colnames(dataset))) {
        pass <- FALSE
        message <- paste(message, "Dataset's 'Gender' column is missed.\n",sep="")
    }    
	
	
    
    Genotype_levels=levels(dataset$Genotype)
    Gender_levels=levels(dataset$Gender) 
    
# Genotype/Gender with at least two data points
    for (i in 1:length(Genotype_levels)){
        GenotypeSubset <- subset(dataset, dataset$Genotype==Genotype_levels[i])
        for (j in 1:length(Gender_levels)){
            nr <- sum(is.finite(GenotypeSubset[GenotypeSubset$Gender==Gender_levels[j],][ , "Gender"]))
            if (nr<2) {
                if (dataset.clean){
                    dataset <- subset(dataset,(dataset$Genotype!=Genotype_levels[i] & dataset$Gender!=Gender_levels[j]))
                    message(paste("Warning: Dataset have been clean: filtered out rows with '", Genotype_levels[i],"'/'",Gender_levels[j],"' combination.\n",sep=""))
                }
                else{
                    pass <- FALSE
                    if (nr==0) nr="no"
                    message_dp <- paste(message_dp,paste("Dataset should have at least two readings for each Genotype/Gender combination. The dataset consists of ", nr, 
														 " reading for '",Genotype_levels[i],"'/'",Gender_levels[j],"' combination.\n",sep=""))
                }
            }    
        }    
        
    }
    
# Renew levels 
    dataset$Gender<-factor(dataset$Gender)
    dataset$Genotype<-factor(dataset$Genotype)
    if ('Batch' %in% colnames(dataset))
	dataset$Batch<-factor(dataset$Batch)
    
    Genotype_levels=levels(dataset$Genotype)
    Gender_levels=levels(dataset$Gender)
    
    genotype_values<-paste(Genotype_levels, collapse="', '" )  
    if (length(Genotype_levels)!=2)  {
        pass <- FALSE
        message <- paste(message,"Dataset's 'Genotype' column has to have two values. At the moment: '",genotype_values,"'\n",sep="")
    }      
    
    if (length(Gender_levels)>2) {
        pass <- FALSE
        message <- paste(message,"Dataset's 'Gender' column has to have one or two values and currently the data has more than two.\n",sep="")
        
    }   
    
    gender_values<-paste(Gender_levels, collapse="', '" )   
    if (!('Female' %in% levels(dataset$Gender)) || !('Male' %in% levels(dataset$Gender))) {
        pass <- FALSE
        message <- paste(message, paste("Dataset's 'Gender' column has '",gender_values,"' values instead of 'Female' and/or 'Male' values. You can specify 'dataset.values.male' and 'dataset.values.female' arguments of the PhenList method to replace values automatically.\n",sep=""),sep="")
    }    
	
    
    message <- paste(message,message_dp,sep="")
    
    
    if (sum(grepl(refGenotype, Genotype_levels, fixed=TRUE))==1)
	dataset$Genotype=relevel(dataset$Genotype, ref=refGenotype)
    else { 
        pass <- FALSE
        message <- paste(message,paste("Dataset with not enough records for statistical analysis with reference genotype '",refGenotype,"'.\n",sep=""))
    }
    if (!pass)
	if (outputMessages){            
		stop(paste("\n",message,sep=""),call. = FALSE) 
	}
	else {
#opt <- options(show.error.messages=FALSE)
#on.exit(options(opt))
		stop()
	}
	
    return(dataset)
    
}