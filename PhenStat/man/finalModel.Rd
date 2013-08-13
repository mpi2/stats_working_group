\name{finalModel}
\alias{finalModel}
\title{Method "finalModel"}
\description{
    Work with PhenTestResult object created by startModel function. 
    The creation of MM final model is based on the significance of different fixed effects, depVariable and equation values 
    stored in phenTestResult object.
}
\usage{
    finalModel(phenList, phenTestResult, outputMessages = TRUE)
}
\arguments{   
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; mandatory argument}
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class; mandatory argument}    
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages; default value TRUE}   

}
\value{
    Return results stored in instance of the \code{\linkS4class{PhenTestResult}} class
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}} and \code{\link{testDataset}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,depVariable="Lean.Mass",equation="withoutWeight",callAll=FALSE)
    result <- finalModel(test,result)
}