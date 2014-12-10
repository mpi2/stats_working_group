\name{parserOutputTFSummary}
\alias{parserOutputTFSummary}
\title{Method "parserOutputTFSummary"}
\description{
    Parsers model output summary from TF framework and returns it in readable named vector format
}
\usage{
    parserOutputTFSummary(phenTestResult)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class that comes from the function \code{\link{testDataset}}; mandatory argument}
}
\value{
    Returns a named vector with model output summary results
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    file <- system.file("extdata", "test6_RR.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Oxr1/Oxr1")
    result <- testDataset(test,
            depVariable="Ca",
            method="TF",
            dataPointsThreshold=2)
    parserOutputTFSummary(result)
}
