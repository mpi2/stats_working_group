# PhenList.R contains PhenList function which construct new PhenList object from components

PhenList <- function(phendata=matrix(0,0,0), refGenotype='+/+', datasetstat=NULL, assay.date.colname=NULL,
        genotype.colname=NULL, gender.colname=NULL, weight.colname=NULL, remove.zeros=FALSE) 

# Construct PhenList object from components with data quality checks

{
    phendata <- phendata[,order(names(phendata))]
    # Rename columns if needed
    if(!is.null(assay.date.colname)) colnames(phendata)[colnames(phendata) == assay.date.colname] <-'Assay.Date'
    if(!is.null(genotype.colname)) colnames(phendata)[colnames(phendata) == genotype.colname] <-'Genotype'
    if(!is.null(gender.colname)) colnames(phendata)[colnames(phendata) == gender.colname] <-'Gender'
    if(!is.null(weight.colname)) colnames(phendata)[colnames(phendata) == weight.colname] <-'Weight'
    
    # Minimum required data
    if (!('Assay.Date' %in% colnames(phendata))) stop("Phenotypic data must have 'Assay.Date' column")
    if (!('Genotype' %in% colnames(phendata))) stop("Phenotypic data must have 'Genotype' column")
    if (!('Gender' %in% colnames(phendata))) stop("Phenotypic data must have 'Gender' column")
    
    # TODO
    # Check and rename when needed genders into Male/Female (M/F, 1/2 etc.)
    
    Genotype_levels=levels(phendata$Genotype)
    
    if(sum(grepl(refGenotype, Genotype_levels, fixed=TRUE))==1){
        phendata$Genotype=relevel(phendata$Genotype, ref=refGenotype)
    }
    else stop(paste("There are not enough records for statistical analysis with reference genotype",refGenotype))
    
    # Statistics
    datasetstat <- data.frame(Variables = colnames(phendata),Numeric = sapply(phendata, is.numeric), 
            NObs = sapply(phendata, function(x) length(na.omit(x))),
            Mean = sapply(phendata, function(x) if(is.numeric(x)) mean(na.omit(x)) else NA),
            StdDev = sapply(phendata, function(x) if(is.numeric(x)) sd(na.omit(x)) else NA),
            Minimum = sapply(phendata, function(x) if(is.numeric(x)) min(na.omit(x)) else NA),
            Maximum = sapply(phendata, function(x) if(is.numeric(x)) max(na.omit(x)) else NA))
    rownames(datasetstat) <- NULL
    
    #    Check variables
    #phendata <- as.matrix(phendata)
    nvar <- ncol(phendata)
    ntags <- nrow(phendata)
    # Column names should be given
    if(nvar>0 && is.null(colnames(phendata))) stop("Phendata must have column names")
    
    if(ntags>0 && is.null(rownames(phendata))) rownames(phendata) <- 1:ntags
    
    x <- new("PhenList",list(phendata=phendata))

    x$datasetstat <- datasetstat
    x$refGenotype <- refGenotype
    x$assay.date.colname <- assay.date.colname
    x$genotype.colname <- genotype.colname
    x$gender.colname <- gender.colname
    x$remove.zeros <- remove.zeros
    
    if(remove.zeros) {
        all.zeros <- rowSums(phendata,na.rm=TRUE)==0
        if(any(all.zeros)) {
            x <- x[!all.zeros,]
            message("Removing ",sum(all.zeros)," rows with all zero counts.")
        }
    }
    
    x
}
