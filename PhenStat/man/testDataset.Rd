\name{testDataset}
\alias{testDataset}
\title{Method "testDataset"}
\description{
    Tests the fixed effects and creates the start model
}
\usage{
    testDataset(phenList, depVariable=NULL, equation="withWeight", outputMessages=TRUE, pThreshold=0.05, method="MM")
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class}       
    \item{depVariable}{a character string defining the dependent variable of interest}
    \item{equation}{a character string defining the equation to use}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages}        
    \item{pThreshold}{a numerical value for the p-value threshold to use in tests}    
    \item{method}{a character string defining the method to use for model building, default value "MM" mixed model}
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenList}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    result <- testDataset(test,depVariable="Bone.Area")
}