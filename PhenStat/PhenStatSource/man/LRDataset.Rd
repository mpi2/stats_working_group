\name{LRDataset}
\alias{LRDataset}
\title{Method "LRDataset"}
\description{  
    Prepares dataset for the LR framework - maps values of dependent variable to 0/1, where 1 is modeled. 
}    
\usage{
    LRDataset(phenList = NULL, depVariable = NULL, abnormalValues =
            c("abnormal", "Abnormal", "TRUE", "deviant"),
            outputMessages = TRUE)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; mandatory argument}       
    \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
    \item{abnormalValues}{a list of abnormal values that will be mapped to 0; mandatory argument}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages; 
        default value TRUE}  
}
\value{
    Returns dataset suitable for LR framework, where provided abonormal values are mapped to 1 and all other values of 
    dependent variable are assumed to be normal and are mapped to 0.
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\examples{
    file <- system.file("extdata", "testLR.csv", package="PhenStat")
    testLR <- PhenList(dataset=read.csv(file),
            testGenotype="Mock")
    testLR2 <- LRDataset(testLR,"V2")
}