\name{testFinalLRModel}
\alias{testFinalLRModel}
\title{Method "testFinalLRModel"}
\description{
    Performs diagnostic tests for Logistic Regression model quality of fit. There are no arguments checks assuming that 
    function is called internally from the "finalLRModel" function. Otherwise should be used with precaution. 
}
\usage{
    testFinalLRModel(phenTestResult)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class that comes 
        from the function \code{testDataset} ; mandatory argument}
}
\value{
    Returns a vector with model fitting results
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    file <- system.file("extdata", "testLR.csv", package="PhenStat")
    testLR <- PhenList(dataset=read.csv(file),
            testGenotype="Mock")
    testLR2 <- LRDataset(testLR,"V2")
    resultLR <- testDataset(testLR2,
            depVariable="V2",
            method="LR")
    testFinalLRModel(resultLR)
}