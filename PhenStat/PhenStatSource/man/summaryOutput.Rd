\name{summaryOutput}
\alias{summaryOutput}
\title{Method "summaryOutput"}
\description{
    Wrapper for the model fitting results. Returns model fitting and testing results in a user friendly format.
}
\usage{
    summaryOutput(phenTestResult, phenotypeThreshold = 0.01)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class that comes from the function \code{\link{testDataset}}; mandatory argument}
    \item{phenotypeThreshold}{a numerical value defining the threshold to use in classificationTag; default value 0.01}
}
\value{
    Returns model fitting results in a text format for the screen
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    # Mixed Model framework
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,
            depVariable="Lean.Mass")
    summaryOutput(result)
    
    # Fishe Exact Test framework
    file <- system.file("extdata", "test_categorical.csv", package="PhenStat")
    test2 <- PhenList(dataset=read.csv(file),
            testGenotype="Aff3/Aff3")
    result2 <- testDataset(test2,
            depVariable="Thoracic.Processes",
            method="FE")  
    summaryOutput(result2)
}