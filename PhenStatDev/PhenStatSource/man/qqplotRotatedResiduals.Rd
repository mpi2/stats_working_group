\name{qqplotRotatedResiduals}
\alias{qqplotRotatedResiduals}
\title{Method "qqplotRotatedResiduals"}
\description{
    Graph function for the Mixed Model framework's results. 
    Creates a Q-Q plot of rotated residuals.
}
\usage{
    qqplotRotatedResiduals(phenTestResult,outputMessages=TRUE)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} 
    class that comes from the function \code{\link{testDataset}}; 
    mandatory argument}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; 
    "TRUE" value to show output messages; default value TRUE}   
}
\value{Creates a plot.}
\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and 
    Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} 
    \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    Houseman, E. A., Ryan, L. M., Coull, B. A. (2004): Cholesky residuals for 
    assessing normal errors in a linear model with correlated outcomes. 
    \emph{Journal of the American Statistical Association} \bold{99}{466}: 
    pg 383-394. Doi 10.1198
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide 
    using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenList}} and \code{\linkS4class{PhenTestResult}}}
\examples{
    file <- system.file("extdata", "test2.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            dataset.colname.weight="Weight.Value",
            testGenotype="Arid4a/Arid4a",
            hemiGenotype="Arid4a/+")
    result <- testDataset(test,
            depVariable="Cl")
    qqplotRotatedResiduals(result)
}