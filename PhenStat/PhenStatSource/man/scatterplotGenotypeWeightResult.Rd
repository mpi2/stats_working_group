\name{scatterplotGenotypeWeightResult}
\alias{scatterplotGenotypeWeightResult}
\title{Method "scatterplotGenotypeWeightResult"}
\description{
    Graph function for the phenotypic dataset. Creates a scatter plot body weight versus dependent variable.
    It is based on PhenTestResult object.
    Both a regression line and a loess line (locally weighted line) is fitted for each genotype.
}
\usage{
    scatterplotGenotypeWeightResult(phenTestResult,graphingName=NULL,outputMessages=TRUE)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class; mandatory argument}
    \item{graphingName}{a character string defining the name to be used on the graph for the dependent variable; mandatory argument}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages; default value TRUE}   
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenList}}}
\examples{
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,
            depVariable="Lean.Mass")              
    scatterplotGenotypeWeightResult(result,
            graphingName="BMC" )
    file <- system.file("extdata", "test4.csv", package="PhenStat")
    test_1sex <- PhenList(dataset=read.csv(file),
            testGenotype="Mysm1/+")
    result_1sex <- testDataset(test_1sex,
            depVariable="Lean.Mass")                 
    scatterplotGenotypeWeightResult(result_1sex,
            graphingName="Lean Mass")
}