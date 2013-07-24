\name{parserOutputSummary}
\alias{parserOutputSummary}
\title{Method "parserOutputSummary"}
\description{
    Parsers model output summary and returns it in readable named vector format
}
\usage{
    parserOutputSummary(phenTestResult)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class}
}
\value{
    Returns a named vector with model output summary results
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
    parserOutputSummary(result)
}
