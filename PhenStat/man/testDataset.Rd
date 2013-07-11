\name{testDataset}
\alias{testDataset}
\title{Tests the fixed effects and creates the start model}
\description{
    Tests the fixed effects and creates the start model
}
\usage{
    testDataset(phenList, interactionMode=TRUE, equation="withWeight", depVariable=NULL, pThreshold=0.05, method="MM")
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class}
    \item{interactionMode}{flag: "FALSE" value for the no interaction messages to the user; "TRUE" value to interact}
    \item{equation}{a character string defining the equation to use}
    \item{depVariable}{a character string defining the depening variable of interest}
    \item{pThreshold}{a numerical value for the p-value threshold to use in tests}    
    \item{method}{a character string defining the method to use for model building, default value "MM" mixed model}
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova}
\seealso{\code{\linkS4class{PhenList}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    result <- testDataset(test,depVariable="Bone.Area")
}