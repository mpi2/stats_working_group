\name{PhenList}
\alias{PhenList}
\title{Method "PhenList"}
\description{
    Function to create data object from the data frame. 
    
    In addition to dependent variable column (the variable of interest) mandatory columns that should present 
    in the data frame are "Genotype" and "Gender". 
    The "Assay.Date" column is used to model "Batch" effect if not specified differently. 
    "Weight" column is used to model body weight effect.
    
    Function creates \code{\linkS4class{PhenList}} object, checks data integrity, renames columns when requested, 
    provides simple statistics about dataset.
    
    Returns an instance of the \code{\linkS4class{PhenList}} object created from the data file. 
    
    Dataset is cleaned to ensure there is a maximum two gender levels and there are exactly two levels for genotype.
    Data cleaning can be switched off by defining the argument "dataset.clean" as FALSE.

}
\usage{
    PhenList(dataset, testGenotype, refGenotype="+/+", hemiGenotype=NULL, 
            outputMessages=TRUE, dataset.clean=TRUE, 
            dataset.colname.batch=NULL, dataset.colname.genotype=NULL, 
            dataset.colname.gender=NULL, dataset.colname.weight=NULL, 
            dataset.values.missingValue=NULL, dataset.values.male=NULL, 
            dataset.values.female=NULL)
    
}
\arguments{
    \item{dataset}{data frame created from file or from another source; mandatory argument}
    \item{testGenotype}{defines the test genotype to be compared to the reference genotype; mandatory argument}
    \item{refGenotype}{defines the reference genotype; assigned default value is "+/+"}    
    \item{hemiGenotype}{defines the genotype value for hemizygous that will be changed to test genotype value}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages; assigned default value is TRUE}
    \item{dataset.clean}{flag: "FALSE" value for no modification of the dataset; "TRUE" value to clean dataset if needed; assigned default value is TRUE}
    \item{dataset.colname.batch}{column name within dataset for the batch effect}
    \item{dataset.colname.genotype}{column name within dataset for the genotype info}
    \item{dataset.colname.gender}{column name within dataset for the gender info}
    \item{dataset.colname.weight}{column name within dataset for the weight info}
    \item{dataset.values.missingValue}{value used as missing value in the dataset}
    \item{dataset.values.male}{value used to label "males" in the dataset}
    \item{dataset.values.female}{value used to label "females" in the dataset}
}
\value{
    Returns an instance of the \code{\linkS4class{PhenList}} class.
}


\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\seealso{\code{\linkS4class{PhenList}}}

\examples{
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    
    file <- system.file("extdata", "test2.csv", package="PhenStat")
    test2 <- PhenList(dataset=read.csv(file),
            testGenotype="Arid4a/Arid4a",
            refGenotype="+/+",
            hemiGenotype="Arid4a/+",
            dataset.colname.weight="Weight.Value")
    
    file <- system.file("extdata", "test3.csv", package="PhenStat")
    test3 <- PhenList(dataset=read.csv(file),
            dataset.clean=TRUE, 
            dataset.values.female=1, 
            dataset.values.male=2, 
            testGenotype="Mysm1/+")
    
}