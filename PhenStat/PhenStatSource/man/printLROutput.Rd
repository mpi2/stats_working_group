\name{printLROutput}
\alias{printLROutput}
\title{Method "printLROutput"}
\description{
    Prints out for linear regression output
}
\usage{    
    printLROutput(phenTestResult,phenotypeThreshold=0.01)
}
\arguments{
      \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class that comes from the function \code{testDataset}; mandatory argument}
    \item{effectValues}{a vector with genotype effect and standard error for the linear regression; mandatory argument}
    \item{phenotypeThreshold}{a numerical value defining the threshold to use in classificationTag; default value 0.01}
}
\value{
    Prints out the  common part of linear regression output.
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
    printLROutput(resultLR)
}
\keyword{internal}