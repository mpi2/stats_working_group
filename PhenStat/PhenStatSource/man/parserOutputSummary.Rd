\name{parserOutputSummary}
\alias{parserOutputSummary}
\title{Method "parserOutputSummary"}
\description{
    Parsers model output summary and returns it in readable named vector format
}
\usage{
    parserOutputSummary(linearRegressionOutput)
}
\arguments{
    \item{linearRegressionOutput}{linear regression output that comes from the MM methods; mandatory argument}
}
\value{
    Returns a named vector with linear regression model output summary results
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,
            depVariable="Lean.Mass")
    linearRegressionOutput <- analysisResults(result)
    parserOutputSummary(linearRegressionOutput)
}
