\name{method}
\alias{method.PhenTestResult}
\title{Retrieve the method used to create the model}
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
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    result <- testDataset(test,depVariable="Bone.Area")
    result <- buildFinalModel(test,result)
    method.PhenTestResult(result)
}

