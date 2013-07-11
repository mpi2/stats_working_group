\name{plotResidualPredicted}
\alias{plotResidualPredicted}
\title{Creates predicted versus residual values plots split by genotype}
\description{
    Creates predicted versus residual values plots split by genotype
}
\usage{
    plotResidualPredicted(phenList, phenTestResult)
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
    plotResidualPredicted(test,result)
}