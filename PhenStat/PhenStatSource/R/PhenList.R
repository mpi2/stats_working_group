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
## PhenList.R contains PhenList and checkDataset functions to construct
## new PhenList object from components and to check dataset integrity
##------------------------------------------------------------------------------
## Construct PhenList object from components with data quality checks
PhenList <- function(dataset, testGenotype, refGenotype='+/+', hemiGenotype=NULL,
        outputMessages=TRUE, dataset.clean=TRUE,
        dataset.colname.batch=NULL, dataset.colname.genotype=NULL,
        dataset.colname.gender=NULL, dataset.colname.weight=NULL,
        dataset.values.missingValue=NULL, dataset.values.male=NULL,
        dataset.values.female=NULL)
{
    dataset <- dataset[,order(names(dataset))]
    
    ## Rename columns if needed
    if (dataset.clean){
        
        if(!is.null(dataset.colname.batch))
        colnames(dataset)[colnames(dataset) == dataset.colname.batch] <-'Batch'
        else {
            if ('Assay.Date' %in% colnames(dataset)){
                colnames(dataset)[colnames(dataset) == 'Assay.Date'] <-'Batch'
                if (outputMessages)
                message("Warning:\nDataset's column 'Assay.Date' has been
                        renamed to 'Batch' and will be used for the batch effect modelling.\n")
            }
            else
            if (length(colnames(dataset)[grep("batch",
                                            tolower(colnames(dataset)))])>0 && outputMessages){
                batch_potential_columns<-paste(colnames(dataset)[grep
                                ("batch", tolower(colnames(dataset)))], collapse="', '" )
                
                message(paste("Warning:\nDataset contains columns that might
                                be used for Batch effect modeling, for instance '",
                                batch_potential_columns,"'.\n",sep=""))
            }
            
        }
        if(!is.null(dataset.colname.genotype))
        colnames(dataset)[colnames(dataset) == dataset.colname.genotype] <-'Genotype'
        if(!is.null(dataset.colname.gender))
        colnames(dataset)[colnames(dataset) == dataset.colname.gender] <-'Gender'
        if(!is.null(dataset.colname.weight))
        colnames(dataset)[colnames(dataset) == dataset.colname.weight] <-'Weight'
        
        ## Replace missing values specified in the user format with NA 
        if(!is.null(dataset.values.missingValue)) 
        dataset[dataset == dataset.values.missingValue] <- NA
        ## Replace empty strings with NA
        dataset[dataset == ""] <- NA
        
        if ('Weight' %in% colnames(dataset)){
            if (is.numeric(dataset$Weight)){
                dataset$Weight<-as.numeric(dataset$Weight)
            }
            else {
                colnames(dataset)[colnames(dataset) == 'Weight'] <-'Weight_labels'
                if (outputMessages)
                message("Warning:\nWeight column values are not numeric. 
                        In order to avoid erroneous execution of statistical 
                        functions column is renamed to 'Weight_labels'.\n")
                
            }
        }
        
        ## Renew levels
        if ('Gender' %in% colnames(dataset))
        dataset$Gender<-factor(dataset$Gender)
        if ('Genotype' %in% colnames(dataset))
        dataset$Genotype<-factor(dataset$Genotype)
        if ('Batch' %in% colnames(dataset))
        dataset$Batch<-factor(dataset$Batch)
        
        # # Replace values for genders with 'Male','Female' if needed
        if(!is.null(dataset.values.female)) 
        levels(dataset$Gender)[levels(dataset$Gender)==dataset.values.female] <- "Female"
        if(!is.null(dataset.values.male)) 
        levels(dataset$Gender)[levels(dataset$Gender)==dataset.values.male] <- "Male"
        
        ## Hemi to test genotype replacement
        if (!is.null(hemiGenotype)) {
            if (length(rownames(dataset[dataset$Genotype==hemiGenotype,]))>0) {
                levels(dataset$Genotype)[levels(dataset$Genotype)==hemiGenotype] <- testGenotype
                if (outputMessages)
                message(paste("Warning:\nHemizygotes '",hemiGenotype,
                                "' have been relabelled to test genotype '",testGenotype,
                                "'.\nIf you don't want this behaviour then don't define 
                                'hemiGenotype' argument.\n",sep=""))
            }
        }
        
        ## Clean genotypes
        if (length(setdiff(rownames(dataset),
                                rownames(dataset[dataset$Genotype %in% c(testGenotype,refGenotype),])))>0){
            dataset <- dataset[dataset$Genotype %in% c(testGenotype,refGenotype),]
            if (outputMessages)
            message(paste("Warning:\nDataset has been cleaned by 
                            filtering out records with genotype value other than test 
                            genotype '", testGenotype,"' or reference genotype '",
                            refGenotype,"'.\n",sep=""))
            
        }
    }
    
    
    ## Clean the empty records  NB - after renaming/cleaning !
    if ('Gender' %in% colnames(dataset))
    dataset<-dataset[dataset$Gender!="",]
    if ('Genotype' %in% colnames(dataset))
    dataset<-dataset[dataset$Genotype!="",]
    if ('Batch' %in% colnames(dataset))
    dataset<-dataset[dataset$Batch!="",]
    
    ## Renew levels
    if ('Gender' %in% colnames(dataset))
    dataset$Gender<-factor(dataset$Gender)
    if ('Genotype' %in% colnames(dataset))
    dataset$Genotype<-factor(dataset$Genotype)
    if ('Batch' %in% colnames(dataset))
    dataset$Batch<-factor(dataset$Batch)
    
    ## CHECKS
    dataset <- checkDataset(dataset, testGenotype, refGenotype, 
            outputMessages, dataset.clean)
    
    
    Genotype_levels <- levels(dataset$Genotype)
    Gender_levels <- levels(dataset$Gender)
    
    ## Calculate statistics
    dataset.stat <- data.frame(
            Variables = colnames(dataset),
            Numeric = sapply(dataset, is.numeric),    
            Continuous = sapply(dataset, function(x) if(is.numeric(x)) {
                        if (length(unique(x))/length(x)>0.05) TRUE else FALSE} else FALSE),
            Levels = sapply(dataset, function(x) length(unique(x)) ),
            NObs = sapply(dataset, function(x) length(na.omit(x))),
            Mean = sapply(dataset, function(x) 
                    if(is.numeric(x)) round(mean(na.omit(x)),digits=2) else NA),
            StdDev = sapply(dataset, function(x) 
                    if(is.numeric(x)) round(sd(na.omit(x)),digits=2) else NA),
            Minimum = sapply(dataset, function(x) 
                    if(is.numeric(x)) round(min(na.omit(x)),digits=2) else NA),
            Maximum = sapply(dataset, function(x) 
                    if(is.numeric(x)) round(max(na.omit(x)),digits=2) else NA))
    
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
#-------------------------------------------------------------------------------
## Check dataset for the minimum required info and additional cleaning steps
checkDataset <- function(dataset, testGenotype, refGenotype="+/+", 
        outputMessages=TRUE, dataset.clean=TRUE)
{
    message <- ""
    message_dp <- ""
    pass <- TRUE
    
    nvar <- ncol(dataset)
    ntags <- nrow(dataset)
    
    ## Column names should be given
    if(nvar>0 && is.null(colnames(dataset))) {
        pass <- FALSE
        message <- paste(message, 
                "\nCheck failed:\nDataset with no column names.\n",sep="")
    }
    
    ## Check for mandatory columns: Genotype and Gender
    if (!('Genotype' %in% colnames(dataset))) {
        pass <- FALSE
        message <- paste(message,
                "\nCheck failed:\nDataset's 'Genotype' column is missed.\n",sep="")
    }
    
    
    if (!('Gender' %in% colnames(dataset))) {
        pass <- FALSE
        message <- paste(message, 
                "\nCheck failed:\nDataset's 'Gender' column is missed.\n",sep="")
    }
    
    ## Check for other columns: Weight and Batch
    if (!('Weight' %in% colnames(dataset)) && outputMessages) {
        message("Warning:\nDataset's 'Weight' column is missed.\n
                You can define 'dataset.colname.weight' argument to specify column 
                for the weight effect modeling. Otherwise you can only use mixed 
                model equation 'withoutWeight'.\n")
    }
    
    if (!('Batch' %in% colnames(dataset))){
        message("Warning:\nDataset's 'Batch' column is missed.\n
                You can define 'dataset.colname.batch' argument to specify column 
                for the batch effect modeling. Otherwise you can only fit a glm.\n")
    }
    
    if (('Gender' %in% colnames(dataset)) && ('Genotype' %in% colnames(dataset))){
        
        Genotype_levels <- levels(dataset$Genotype)
        Gender_levels <- levels(dataset$Gender)
        
        ## Genotype/Gender combinations with less than two data points
        combinations_list <- "" # String with Genotype/Gender (count) in initial dataset
        
        for (i in 1:length(Genotype_levels)){
            GenotypeSubset <- subset(dataset, dataset$Genotype==Genotype_levels[i])
            for (j in 1:length(Gender_levels)){
                nr <- sum(is.finite(GenotypeSubset[GenotypeSubset$Gender==Gender_levels[j],][ , "Gender"]))
                combinations_list<-paste(combinations_list,paste("'", 
                                Genotype_levels[i],"'/'",Gender_levels[j],
                                "' (",nr,"), ",sep=""),sep="")
            }
            
        }
        
        ## String with Genotype/Gender (count), where Gender would be filtered out
        filtered_list_combinations <- ""
        dataset_filtered <- dataset
        for (i in 1:length(Genotype_levels)){
            GenotypeSubset <- subset(dataset, dataset$Genotype==Genotype_levels[i])
            for (j in 1:length(Gender_levels)){
                nr <- sum(is.finite(GenotypeSubset[GenotypeSubset$Gender==Gender_levels[j],][ , "Gender"]))
                
                ## There are combinations with less than two data points
                if (nr<2) {
                    filtered_list_combinations <- paste(
                            filtered_list_combinations,
                            paste("'", Genotype_levels[i],"'/'",Gender_levels[j],
                                    "' (",nr,"), ",sep=""),sep="")
                    
                    if (dataset.clean){
                        ## If you have data in one genotype for both genders 
                        ## but not in the other then you have to revert to a 
                        ## one gender analysis.
                        subset_to_filter <- subset(dataset_filtered,
                                (dataset_filtered$Gender==Gender_levels[j]))
                        
                        dataset_filtered <- dataset_filtered[setdiff(
                                        rownames(dataset_filtered),rownames(subset_to_filter)),]
                        
                    }
                    else{
                        pass <- FALSE
                    }
                }
            }
            
        }
        
        
        if (nchar(filtered_list_combinations)>2){
            combinations_list <- substr(combinations_list, 1, 
                    nchar(combinations_list)-2)
            
            filtered_list_combinations <- substr(filtered_list_combinations, 
                    1, nchar(filtered_list_combinations)-2)
            
            if (dataset.clean){
                if (outputMessages)
                message(paste("Warning:\nSince dataset has to have at least 
                                two data points for each genotype/gender combination and 
                                there are not enough records for the combination(s): ",
                                filtered_list_combinations,", appropriate gender records 
                                have been filtered out from the dataset.\n",sep=""))
            }
            else
            message_dp <- paste("\nCheck failed:\nDataset should have at 
                    least two data points for each genotype/gender combination. 
                    At the moment there are no enough data points for the following 
                    combination(s): ",filtered_list_combinations,".\n",sep="")
        }
        
        dataset <- dataset_filtered
        
        ## Renew levels
        dataset$Gender <- factor(dataset$Gender)
        dataset$Genotype <- factor(dataset$Genotype)
        
        if ('Batch' %in% colnames(dataset))
        dataset$Batch <- factor(dataset$Batch)
        
        
        Genotype_levels <- levels(dataset$Genotype)
        Gender_levels <- levels(dataset$Gender)
        
        ## INFO about genotype and gender levels
        genotype_values <- paste(Genotype_levels, collapse="', '" )
        gender_values <- paste(Gender_levels, collapse="', '" )
        if (outputMessages){
            message("Information:\nDataset's 'Genotype' 
                    column has following values: '",genotype_values,"'\n",sep="")
            message("Information:\nDataset's 'Gender' 
                    column has following value(s): '",gender_values,"'\n",sep="")
        }
        
        ## Check of genotype and gender levels after cleaning
        if (length(Genotype_levels)!=2)  {
            pass <- FALSE
            message <- paste(message,"\nCheck failed:\nDataset's 'Genotype' 
                    column has to have two values.\nYou can define 'testGenotype' and
                    'refGenotype' arguments to automatically filter out records with 
                    genotype values other than specified. Alternatively you can define 
                    'hemiGenotype' and 'testGenotype' arguments to relabel hemizygotes 
                    to homozygotes.\n",sep="")
        }
        
        
        if (!(length(Gender_levels) %in% c(1,2))) {
            pass <- FALSE
            message <- paste(message,"\nCheck failed:\nDataset's 'Gender' 
                    column has to have one or two values and currently the data has 
                    more than two.\n",sep="")
            
        }
        
        ## Check for gender levels - we want to have 'Female' and/or 'Male' only
        wrong_gender_levels <- setdiff(Gender_levels,c("Female","Male"))
        wrong_gender_values<-paste(wrong_gender_levels, collapse="', '" )
        
        if (!length(wrong_gender_levels)==0){
            pass <- FALSE
            if (length(Gender_levels)<=2)
            message <- paste(message, paste("\nCheck failed:\nDataset's 
                            'Gender' column has '",gender_values,"' values instead of 
                            'Female' and/or 'Male' values. You can define 
                            'dataset.values.male' and 'dataset.values.female' 
                            arguments to replace those values automatically.\n",sep=""),
                    sep="")
            else
            message <- paste(message, paste("\nCheck failed:\nDataset's 
                            'Gender' column has '",gender_values,"' values instead of 
                            'Female' and/or 'Male' values only. 
                            Please delete records with gender(s) '",wrong_gender_values,"' 
                            from the dataset.\n",sep=""),sep="")
                }
                
                
                message <- paste(message,message_dp,sep="")
                
                ## Check for reference genotype records
                if (sum(grepl(refGenotype, Genotype_levels, fixed=TRUE))==1)
                dataset$Genotype=relevel(dataset$Genotype, ref=refGenotype)
                else {
                    pass <- FALSE
                    message <- paste(message,paste("\nCheck failed:\nDataset with not 
                                    enough records for statistical analysis with reference genotype '",
                                    refGenotype,"'.\n",sep=""))
                }
                
                ## Check for test genotype records
                if (!(sum(grepl(testGenotype, Genotype_levels, fixed=TRUE))==1)){
                    pass <- FALSE
                    message <- paste(message,paste("\nCheck failed:\nDataset 
                                    with not enough records for statistical analysis with test 
                                    genotype '",testGenotype,"'.\n",sep=""))
                }
            }
            
            if (!pass){
                if (outputMessages){
                    message(paste("********* Errors start *********\n",message,sep=""))
                    message("********* Errors end ***********")
                }
                opt <- options(show.error.messages=FALSE)
                on.exit(options(opt))
                stop()
            }
            
            return(dataset)
}
##------------------------------------------------------------------------------