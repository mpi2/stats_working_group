\name{buildStartModel}
\alias{buildStartModel}
\title{Method "buildStartModel"}
\description{
    If someone would like to assign other TRUE/FALSE values to effects of the model then start model 
    is built by using this function. There are no dataset checks or arguments checks  as it is assumed that buildStartModel
    function is called internally from the buildFinalModel function. Therefore this function should be used with precaution. 
}
\usage{    
    buildStartModel(object, equation, depVariable, keepList)
}
\arguments{
    \item{object}{instance of the \code{\linkS4class{PhenList}} class}
    \item{equation}{a character string defining the equation to use}
    \item{depVariable}{a character string defining the dependent variable of interest}
    \item{keepList}{a logical vector defining the significance of different effects: keep_batch, keep_equalvar,
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
    result <- buildStartModel(test,equation="withWeight",depVariable="Bone.Area",keepList=c(TRUE,TRUE,TRUE,TRUE,FALSE))
}