\name{JSONOutput}
\alias{JSONOutput}
\title{Method "JSONOutput"}
\description{
    Wrapper for the model fitting results. Returns model fitting and testing results in a JSON format. 
    Assumes that modelling results are stored in the \code{\linkS4class{PhenTestResult}} object 
    (output from function \code{\link{testDataset}}).
}
\usage{
    JSONOutput(phenTestResult, phenotypeThreshold = 0.01)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class that comes from the function \code{\link{testDataset}}; mandatory argument}
    \item{phenotypeThreshold}{a numerical value defining the threshold to use in classificationTag; default value 0.01}
}
\value{
    Returns a vector with model fitting results in JSON format: model output summary and other values
}

\author{Natalja Kurbatova}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    # Mixed Model framework
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,
            depVariable="Lean.Mass")
    JSONOutput(result)
    
    # Fishe Exact Test framework
    file <- system.file("extdata", "test_categorical.csv", package="PhenStat")
    test2 <- PhenList(dataset=read.csv(file),
            testGenotype="Aff3/Aff3")
    result2 <- testDataset(test2,
            depVariable="Thoracic.Processes",
            method="FE")  
    JSONOutput(result2)
}