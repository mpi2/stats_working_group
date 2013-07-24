\name{qqplotRotatedResiduals}
\alias{qqplotRotatedResiduals}
\title{Method "qqplotRotatedResiduals"}
\description{
    Creates a Q-Q plot of rotated residuals
}
\usage{
    qqplotRotatedResiduals(phenList, phenTestResult, keep_batch)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class}
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class}
    \item{keep_batch}{a logical value defining the significance of batch effect}
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenList}} and \code{\linkS4class{PhenTestResult}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MAAJ_CCdata.csv"),dataset.colname.weight="Weight.Value",
            dataset.clean=TRUE,testGenotype="Arid4a/Arid4a",refGenotype="+/+",hemiGenotype="Arid4a/+")
    result <- testDataset(test,depVariable="Cl")
    result <- buildFinalModel(test,result)
    qqplotRotatedResiduals(test,result)
}