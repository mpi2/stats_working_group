\name{buildFinalModel}
\alias{buildFinalModel}
\title{Method "buildFinalModel"}
\description{
    Build final mixed model
}
\usage{
    buildFinalModel(phenList, phenTestResult=NULL, depVariable=NULL, equation=NULL, outputMessages=TRUE, keepList=NULL)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class}
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class}    
    \item{depVariable}{a character string defining the dependent variable of interest}
    \item{equation}{a character string defining the equation to use}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages}   
    \item{keepList}{a logical vector defining the significanse of different effects: keep_batch, keep_equalvar,
        keep_weight, keep_gender, keep_interaction}
}
\value{
    Returns results stored in instance of the \code{\linkS4class{PhenTestResult}} class
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
    buildFinalModel(test,result)
}