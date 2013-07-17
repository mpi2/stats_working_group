\name{testFinalModel}
\alias{testFinalModel}
\title{Method "testFinalModel"}
\description{
    Performs diagnostic tests for mixed model quality of fit
}
\usage{
    testFinalModel(object,result)
}
\arguments{
    \item{object}{instance of the \code{\linkS4class{PhenList}} class}
    \item{result}{instance of the \code{\linkS4class{PhenTestResult}} class}
}
\value{
    Returns a vector with model fitting results
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
    testFinalModel(test,result)
}