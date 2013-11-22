\name{method}
\alias{method.PhenTestResult}
\title{Method "method"}
\description{
    Retrieve the method used to create the model: gls, lme.
}
\usage{
    \method{method}{PhenTestResult}(x)
}
\arguments{
    \item{x}{an object of class \code{PhenTestResult}}
}
\value{
    Character describing the method used
}
\author{Natalja Kurbatova}
\examples{
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,
            depVariable="Lean.Mass")
    method.PhenTestResult(result)
}

\keyword{internal}