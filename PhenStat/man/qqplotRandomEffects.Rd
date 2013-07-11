\name{qqplotRandomEffects}
\alias{qqplotRandomEffects}
\title{Creates a Q-Q plot of blups}
\description{
    Creates a Q-Q plot of blups
}
\usage{
    qqplotRandomEffects(phenList, phenTestResult, keep_batch)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class}
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class}
    \item{keep_batch}{a logical value defining the significance of batch effect}
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
    qqplotRandomEffects(test,result)
}