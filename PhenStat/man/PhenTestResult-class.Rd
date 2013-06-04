\name{PhenTestResult}
\alias{PhenTestResult}
\title{Returns an instance of the PhenTestResult object.}
\description{
    Creates PhenTestResult object, based on test results or model building results. 
}
\usage{
    PhenTestResult(modelOutput=NULL, depVariable=NULL, equation="withWeight", batchEffect=NULL, varianceEffect=NULL, 
                    genderEffect=NULL, interactionEffect=NULL, interactionTest=NULL, weightEffect=NULL, outputLength=NULL,
                    pThreshold=0.05, model.formula=NULL, model.null=NULL, 
                    model.genotype=NULL, model=NULL, genotypeEffect=NULL, MM_fitquality=NULL)
}
\arguments{
    \item{dataFrame}{data frame created from file or from another source}
}
\value{
    Returns an instance of the \code{\linkS4class{PhenTestResult}} class.
}
\references{
    Publication placeholder

}
\author{Natalja Kurbatova}

\seealso{\code{\linkS4class{PhenTestResult}} and \code{\linkS4class{PhenList}}}

\examples{
    test <- PhenList(phendata=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    result <- testDataset(test,equation="withoutWeight",depVariable="Bone.Area")
}