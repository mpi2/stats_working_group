\name{finalLRModel}
\alias{finalLRModel}
\title{Method "finalLRModel"}
\description{
    This is an internal function run within LR framework. 
    It completes the final stage of the LR framework, which builds the final model and estimates effects.  
    As an internal function, it doesn't include extensive error testing of inputs. Please use cautiously if calling directly.  
    
    Works with \code{\linkS4class{PhenTestResult}} object created by \code{\link{startLRModel}} function. 
    The creation of LR final model is based on the significance of different fixed effects, depVariable and equation values 
    stored in \code{\linkS4class{PhenTestResult}} object.
}
\usage{
    finalLRModel(phenTestResult, outputMessages = TRUE)
}
\arguments{   
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} 
    class that comes from the function \code{\link{testDataset}}; mandatory argument}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages; default value TRUE}   
    
}
\value{
    Returns results stored in instance of the \code{\linkS4class{PhenTestResult}} class
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}} and \code{\link{testDataset}}}
\examples{
    file <- system.file("extdata", "testLR.csv", package="PhenStat")
    testLR <- PhenList(dataset=read.csv(file),
            testGenotype="Mock")
    testLR2 <- LRDataset(testLR,"V2")
    
    # when "testDataset" function's argument "callAll" is set to FALSE 
    # only "startLRModel" function is called - the first step of LR framework
    resultLR <- testDataset(testLR2,
            depVariable="V2",
            callAll=FALSE,
            method="LR")
    # print out formula that has been created
    analysisResults(resultLR)$model.formula.genotype
    # print out batch effect's significance 
    analysisResults(resultLR)$model.effect.batch
    resultLR <- finalLRModel(resultLR)
}