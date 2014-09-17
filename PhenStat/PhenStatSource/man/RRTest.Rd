\name{RRTest}
\alias{RRTest}
\title{Method "RRTest"}
\description{
    This is an internal function run within RR framework. It performs Reference 
    Ranges Plus test and after that 
    Fisher Exact test on calculated count matrices. 
    As an internal function, it doesn't include extensive error testing of 
    inputs. 
    Please use cautiously if calling directly.  
    
    Works with \code{\linkS4class{PhenList}} object created by 
    \code{\link{PhenList}} function.     
}
\usage{
    RRTest(phenList, depVariable, outputMessages=TRUE, naturalVariation=95, 
    controlPointsThreshold=60)
}
\arguments{   
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; 
    mandatory argument}
    \item{depVariable}{a character string defining the dependent variable of 
    interest; mandatory argument}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; 
    "TRUE" value to show output messages; default value TRUE}  
    \item{naturalVariation}{number defining the natural variation range in
    percents for normal values; default value 95}
    \item{controlPointsThreshold}{number defining how many control points is
    needed to perform RR plus method; default value 60} 
    
}
\value{Returns results stored in instance of the 
    \code{\linkS4class{PhenTestResult}} class
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and 
    Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} 
    \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide 
    using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenList}} and \code{\link{testDataset}}}
\examples{
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    # "RRTest" function is called from "testDataset" function
    result <- testDataset(test,
            depVariable="Lean.Mass",
            method="RR")
    # RR test results can be printed out using function "summaryOutput"
    summaryOutput(result)
}