\name{summaryOutput}
\alias{summaryOutput}
\title{Method "summaryOutput"}
\description{
    Wrapper for the model fitting results. Returns model fitting and testing results in a user friendly format
}
\usage{
    summaryOutput(phenTestResult, phenotypeThreshold = 0.01)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class}
    \item{phenotypeThreshold}{a numerical value defining the threshold to use in classification tag}
}
\value{
    Returns model fitting results in a text format
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
    summaryOutput(result)
}