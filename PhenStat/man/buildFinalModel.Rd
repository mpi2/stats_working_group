\name{buildFinalModel}
\alias{buildFinalModel}
\title{Build final mixed model}
\description{
    Build final mixed model
}
\usage{
    buildFinalModel(phenList, phenTestResult=NULL, equation=NULL, depVariable=NULL, pThreshold=0.05, keepList=NULL)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class}
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class}
    \item{equation}{a character string defining the equation to use}
    \item{depVariable}{a character string defining the depening variable of interest}
    \item{pThreshold}{a numerical value for the p-value threshold to use in tests}    
    \item{keepList}{a logical vector defining the significanse of different effects: keep_batch, keep_equalvar,
        keep_weight, keep_gender, keep_interaction}
}
\value{
    Returns results stored in instance of the \code{\linkS4class{PhenTestResult}} class
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    result <- testDataset(test,depVariable="Bone.Area")
    result <- buildFinalModel(test,result)
    buildFinalModel(test,result)
}