\name{startLRModel}
\alias{startLRModel}
\title{Method "startLRModel"}
\description{      
    This is an internal function run within LR framework. It completes the testing stage of which effects are significant.  
    As an internal function, it doesn't include extensive error testing of inputs. Please use cautiously if calling directly.  
    
    It creates start model and modify it after testing of different hypothesis.
    
    The tested effects are: 
    
    -batch effect (TRUE if batch variation is significant, FALSE if not), though it is never used in final LR model; 
    
    -interaction effect (TRUE if genotype by sex interaction is significant), 
    
    -sex effect (TRUE if sex is significant), 
    
    -weight effect and variance effect are not tested and used in LR model. 
    
}    
\usage{
    startLRModel(phenList, depVariable, outputMessages=TRUE, pThreshold=0.05)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; mandatory argument}       
    \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages
        ; default value TRUE}        
    \item{pThreshold}{a numerical value for the p-value threshold used to determine which fixed effects to keep 
        in the model, default value 0.05}    
}
\value{
    Returns results stored in instance of the \code{\linkS4class{PhenTestResult}} class
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenList}}}
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
}