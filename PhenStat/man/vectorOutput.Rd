\name{vectorOutput}
\alias{vectorOutput}
\title{Method "vectorOutput"}
\description{
    Wrapper for the model fitting results. Return model fitting and testing results in a vector format. 
    Assume that modelling results are stored in the phenTestResult object 
    (output from functions testDataset and buildFinalModel).
}
\usage{
    vectorOutput(phenTestResult)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class}
}
\value{
    Return a vector with model fitting results: model output summary and other values
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,depVariable="Lean.Mass")
    vectorOutput(result)
}