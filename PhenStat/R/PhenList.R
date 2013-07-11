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


# PhenList.R contains PhenList and checkDataset functions to construct new PhenList object from components and to check 
# dataset integrity

PhenList <- function(dataset=matrix(0,0,0), interactionMode=TRUE, refGenotype='+/+', dataset.stat=NULL, dataset.colname.batch=NULL,
        dataset.colname.genotype=NULL, dataset.colname.gender=NULL, dataset.colname.weight=NULL, dataset.values.missingValue=NULL, 
        dataset.values.male=NULL, dataset.values.female=NULL, dataset.clean=FALSE, testGenotype=NULL, hemGenotype=NULL) 

# Construct PhenList object from components with data quality checks

# TODO 
# testGenotype
# clean all other records (with other genotypes)

{
    dataset <- dataset[,order(names(dataset))]
    
    
    # Rename columns if needed    
    if(!is.null(dataset.colname.batch)) colnames(dataset)[colnames(dataset) == dataset.colname.batch] <-'Batch'
    else  colnames(dataset)[colnames(dataset) == 'Assay.Date'] <-'Batch'
    if(!is.null(dataset.colname.genotype)) colnames(dataset)[colnames(dataset) == dataset.colname.genotype] <-'Genotype'
    if(!is.null(dataset.colname.gender)) colnames(dataset)[colnames(dataset) == dataset.colname.gender] <-'Gender'
    if(!is.null(dataset.colname.weight)) colnames(dataset)[colnames(dataset) == dataset.colname.weight] <-'Weight'
    
    # Replace missing values specified in the user format with NA if needed 
    if(!is.null(dataset.values.missingValue)) dataset[dataset == dataset.values.missingValue] <- NA 
    
    # Replace values for genders with 'Male','Female' if needed 
    if(!is.null(dataset.values.female)) levels(dataset$Gender)[levels(dataset$Gender)==dataset.values.female] <- "Female"
    if(!is.null(dataset.values.male)) levels(dataset$Gender)[levels(dataset$Gender)==dataset.values.male] <- "Male"
    

    if (dataset.clean && !is.null(hemGenotype) && !is.null(testGenotype)) { 
        levels(dataset$Genotype)[levels(dataset$Genotype)==hemGenotype] <- testGenotype
        if (interactionMode)
            message(paste("Warning: Hemizygotes '",hemGenotype,"' have been relabled to homozygotes '",testGenotype,"'. If you don't want this behaviour please delete the hemizygotes records from the dataset.",sep=""))
    }
    if (dataset.clean && !is.null(testGenotype) && length(levels(dataset$Genotype))>2) {
        dataset <- dataset[-(dataset$Genotype!=testGenotype || dataset$Genotype!=refGenotype),]
        if (interactionMode)
            message(paste("Warning: Dataset has been cleaned, filtered out rows with genotype value other than ", testGenotype,"or",refGenotype))
    }
     
    # Clean the empty records   
    dataset<-dataset[dataset$Gender!="",]
    dataset<-dataset[dataset$Genotype!="",]
    dataset<-dataset[dataset$Batch!="",]
    
    # Renew levels 
    dataset$Gender<-factor(dataset$Gender)
    dataset$Genotype<-factor(dataset$Genotype)
    dataset$Batch<-factor(dataset$Batch)

    dataset <- checkDataset(dataset, interactionMode, refGenotype, dataset.clean)
    
    Genotype_levels=levels(dataset$Genotype)
    Gender_levels=levels(dataset$Gender)
    
    #TODO Reset the levels!
    
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
    x$hemGenotype <- hemGenotype
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

checkDataset <- function(dataset, interactionMode=TRUE, refGenotype="+/+", dataset.clean=FALSE)

# Check dataset for the minimum required info

{
    message <- " "
    pass <- TRUE
    
    
    nvar <- ncol(dataset)
    ntags <- nrow(dataset)

    # Column names should be given
    if(nvar>0 && is.null(colnames(dataset))) {
        pass <- FALSE
        message <- paste(message, "Dataset with no column names")
    }    
    
    if(ntags>0 && is.null(rownames(dataset))) rownames(dataset) <- 1:ntags
    
    # Minimum required data
    #if (!('Batch' %in% colnames(dataset))){
    #    pass <- FALSE
    #    message <- paste(message, "Dataset's 'Batch' column is missed\n")
    #}    
    
    if (!('Genotype' %in% colnames(dataset))) {
        pass <- FALSE
        message <- paste(message,"Dataset's 'Genotype' column is missed\n")
    }
    
    # What about datasets without Gender info???
    if (!('Gender' %in% colnames(dataset))) {
        pass <- FALSE
        message <- paste(message, "Dataset's 'Gender' column is missed\n")
    }    
    
    Genotype_levels=levels(dataset$Genotype)
    Gender_levels=levels(dataset$Gender) 
    # Genotype/Gender with at least two data points
    # Genotype/Gender with at least two data points
    for (i in 1:length(Genotype_levels)){
        GenotypeSubset <- subset(dataset, dataset$Genotype==Genotype_levels[i])
        for (j in 1:length(Gender_levels)){
            nr <- sum(is.finite(GenotypeSubset[GenotypeSubset$Gender==Gender_levels[j],][ , "Gender"]))
            #message(paste(Genotype_levels[i],Gender_levels[j],nr))
            if (nr<2) {
                if (dataset.clean){
                    dataset <- subset(dataset,(dataset$Genotype!=Genotype_levels[i] & dataset$Gender!=Gender_levels[j]))
                    message(paste("Warning: Dataset have been clean: filterd out rows with ", Genotype_levels[i],"/",Gender_levels[j],"combination"))
                }
                else{
                    pass <- FALSE
                    if (nr==0) nr="no"
                    message <- paste(message,paste("Dataset should have at least two readings for each Genotype/Gender combination. The dataset consists of", nr, 
                                    "reading for",Genotype_levels[i],"/",Gender_levels[j],"combination.\n"))
                }
            }    
        }    
        
    }
    
    # Renew levels 
    dataset$Gender<-factor(dataset$Gender)
    dataset$Genotype<-factor(dataset$Genotype)
    dataset$Batch<-factor(dataset$Batch)
    
    Genotype_levels=levels(dataset$Genotype)
    Gender_levels=levels(dataset$Gender)
    
    if (length(Genotype_levels)!=2)  {
        pass <- FALSE
        message <- paste(message,"Dataset's 'Genotype' column have to have two values\n")
    }      
    
    if (length(Gender_levels)>2) {
        pass <- FALSE
        message <- paste(message,"Dataset's 'Gender' column have to have one or two values\n")
    }    
    
    
    
    
    if (sum(grepl(refGenotype, Genotype_levels, fixed=TRUE))==1)
    dataset$Genotype=relevel(dataset$Genotype, ref=refGenotype)
    else { 
        pass <- FALSE
        message <- paste(message,paste("Dataset with not enough records for statistical analysis with reference genotype",refGenotype))
    }
    if (!pass) 
        if (interactionMode)
            stop(paste("Error(s):", message)) 
        else {
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt))
            stop()
        }

    return(dataset)
    
}
