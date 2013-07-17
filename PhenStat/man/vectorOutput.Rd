\name{vectorOutput}
\alias{vectorOutput}
\title{Method "vectorOutput"}
\description{
    Wrapper for the model fitting results. Returns model fitting and testing results in a vector format 
}
\usage{
    vectorOutput(phenTestResult)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class}
}
\value{
    Returns a vector with model fitting results: model output summary and other values
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    result <- testDataset(test,depVariable="Bone.Area")
    result <- buildFinalModel(test,result)
    vectorOutput(result)
}