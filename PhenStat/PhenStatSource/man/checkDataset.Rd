\name{checkDataset}
\alias{checkDataset}
\title{Method "checkDataset"}
\description{
    Checks dataset for the minimum required info:
    
    1. Column names should present
    
    2. Genotype column should present
    
    3. Gender column should present
    
    4. Two data points for each Genotype/Gender combination
    
    5. Number of Genotype levels should be 2
    
    6. Number of Gender levels should be 1 or 2
    
    7. Gender levels have to be "Female" and/or "Male"
    
    8. Records with reference genotype should be in the dataset
    
    9. Records with test genotype should be in the dataset
    
    Perform the following additional checks: 
    
    - presence of Weight column, 
    
    - presence of Batch column (Assay.Date).
    
    Warning given in case of missed data indicating that you can only fit a glm or to use MM equation "withoutWeight".
    
    Function \code{\link{checkDataset}} is called from \code{\link{PhenList}} function. 
}
\usage{
    checkDataset(dataset, testGenotype, refGenotype="+/+", 
            outputMessages=TRUE, dataset.clean=TRUE) 
    
}
\arguments{
    \item{dataset}{data frame created from file or from another source; mandatory argument}
    \item{testGenotype}{defines the test genotype to be compared to the reference genotype; mandatory argument}
    \item{refGenotype}{defines the reference genotype; assigned default value is "+/+"}    
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages; assigned default value is TRUE}
    \item{dataset.clean}{flag: "FALSE" value for no modification of the dataset; "TRUE" value to clean dataset if needed; assigned default value is TRUE}
      
}
\value{
    Returns an instance of the \code{\linkS4class{PhenList}} class.
}


\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\seealso{\code{\linkS4class{PhenTestResult}}}

\examples{    
    # "checkDataset" is called internally from "PhenList" function
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")    
    
}