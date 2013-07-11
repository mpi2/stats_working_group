\name{buildStartModel}
\alias{buildStartModel}
\title{Build start model}
\description{
    If someone would like to assign other TRUE/FALSE values to effects of the model then start model 
    is built by using this function 
}
\usage{    
    buildStartModel(object, equation = "withWeight", depVariable, pThreshold=0.05, keepList)
}
\arguments{
    \item{object}{instance of the \code{\linkS4class{PhenList}} class}
    \item{equation}{a character string defining the equation to use}
    \item{depVariable}{a character string defining the depening variable of interest}
    \item{pThreshold}{a numerical value for the p-value threshold to use in tests}    
    \item{keepList}{a logical vector defining the significanse of different effects: keep_batch, keep_equalvar,
        keep_weight, keep_gender, keep_interaction}
}
\value{
    Returns start model
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    model <- buildStartModel(test,equation="withWeight",depVariable="Bone.Area",keepList=c(TRUE,TRUE,TRUE,TRUE,FALSE))
}