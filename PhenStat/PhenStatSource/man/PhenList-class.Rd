\name{PhenList-class}
\docType{class}
\alias{PhenList-class}

\title{Class "PhenList"}

\description{A list-based S4 class for storing phenotypic data. 
Helps to support data integrity checks and statistics calculation.
The PhenList object can be created by using function \code{\link{PhenList}}.}

\section{Explore PhenList object}{
    PhenList object instance contains the following sections:
    
    1. "dataset" where cleaned and checked dataset is stored: 
    \code{phenList$dataset}
    
    2. "testGenotype" where the genotype level to test is stored: 
    \code{phenList$testGenotype}
    
    3. "refGenotype" where reference genotype  value is stored with 
    default value set to "+/+":  
    
    \code{phenList$refGenotype}
    
    4. "dataset.stat" where simple dataset statistics is stored: 
    \code{phenList$dataset.stat}
    
    5. "hemiGenotype" if defined contains the genotype value for hemizygous: 
    \code{phenList$hemiGenotype}
    
    6. "dataset.clean" flag value is stored,  see \code{\link{PhenList}} 
    for more details:
    
    \code{phenList$dataset.clean}
    
    7. "dataset.colname" if defined contains dataset column names that have 
    been renamed:
    
        - \code{phenList$dataset.colname.batch} column name for batch values
       
        - \code{phenList$dataset.colname.genotype} column name for genotype 
        values
    
        - \code{phenList$dataset.colname.sex} column name for sex values
    
        - \code{phenList$dataset.colname.weight} column name for weight values
        
    8. "dataset.values" if defined contains dataset values that have been 
    changed during dataset cleaning process:
    
        - \code{phenList$dataset.values.missingValue} value used as missing 
        value in the orginal dataset
        
        - \code{phenList$dataset.values.male} value used to label "males" in 
        the original dataset
        
        - \code{phenList$dataset.values.female} value used to label "females" 
        in the original dataset
    

}
\value{PhenList class}
\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and 
    Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} 
    \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide 
    using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}

\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\seealso{\code{\link{PhenList}}}

\examples{
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    class(test)
    
    file <- system.file("extdata", "test2.csv", package="PhenStat")
    test2 <- PhenList(dataset=read.csv(file),
            testGenotype="Arid4a/Arid4a",
            refGenotype="+/+",
            hemiGenotype="Arid4a/+",
            dataset.colname.weight="Weight.Value")
    test2$dataset.stat
    test2$testGenotype
    test2$refGenotype
    test2$hemiGenotype
    test2$dataset.colname.weight
}


\keyword{classes}