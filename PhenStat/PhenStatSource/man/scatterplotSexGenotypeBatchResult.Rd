\name{scatterplotSexGenotypeBatchResult}
\alias{scatterplotSexGenotypeBatchResult}
\title{Method "scatterplotSexGenotypeBatchResult"}
\description{
    Graph function for the phenotypic dataset. Creates a scatterplot split by sex, genotype and batch. refGenotype data points are shown in black and the testGenotype data points are shown in red.
    It is based on PhenTestResult object.
    
    Note: the batches are not ordered with time but allow assessment of how the testGenotype data lie relative to the variation within the refGenotype data.
    
}
\usage{
    scatterplotSexGenotypeBatchResult(phenTestResult,graphingName=NULL,outputMessages=TRUE)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class; mandatory argument}
    \item{graphingName}{a character string defining the name to be used on the graph for the dependent variable; mandatory argument}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages; default value TRUE}   
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
    result <- testDataset(test,
            depVariable="Lean.Mass")        
    # box plot for dataset with two sexes: males and females
    scatterplotSexGenotypeBatchResult(result,
            graphingName="BMC" )
    
    file <- system.file("extdata", "test4.csv", package="PhenStat")
    test_1sex <- PhenList(dataset=read.csv(file),
            testGenotype="Mysm1/+")
    result_1sex <- testDataset(test_1sex,
            depVariable="Lean.Mass")        
    # box plot for females only dataset
    scatterplotSexGenotypeBatchResult(result_1sex, 
            graphingName="Lean Mass")
}