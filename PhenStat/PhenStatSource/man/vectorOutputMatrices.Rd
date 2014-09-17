\name{vectorOutputMatrices}
\alias{vectorOutputMatrices}
\title{Method "vectorOutputMatrices"}
\description{
    Returns count matrices from Fisher Exact Test framework in a vector form.
    
    Note: Maximal supported number of dependent variable levels is 10.
}
\usage{
    vectorOutputMatrices(phenTestResult,outputMessages=TRUE)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} 
    class that comes from the function \code{\link{testDataset}}; 
    mandatory argument}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; 
    "TRUE" value to show output messages; default value TRUE}   
}
\value{Returns a vector with count values.
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and 
    Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} 
    \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide 
    using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    file <- system.file("extdata", "test_categorical.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Aff3/Aff3")
    result <- testDataset(test,
            depVariable="Thoracic.Processes",
            method="FE")  
    vectorOutputMatrices(result)
}