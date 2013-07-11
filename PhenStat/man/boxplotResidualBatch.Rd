\name{boxplotResidualBatch}
\alias{boxplotResidualBatch}
\title{Creates a box plot with residue versus batch split by genotype}
\description{
    Creates a box plot with residue versus batch split by genotype
}
\usage{
    boxplotResidualBatch(phenList, phenTestResult)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class}
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class}
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova}
\seealso{\code{\linkS4class{PhenList}} and \code{\linkS4class{PhenTestResult}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    result <- testDataset(test,depVariable="Bone.Area")
    result <- buildFinalModel(test,result)
    boxplotResidualBatch(test,result)
}