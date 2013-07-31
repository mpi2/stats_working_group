\name{buildFisherExactTest}
\alias{buildFisherExactTest}
\title{Method "buildFisherExactTest"}
\description{
    Create two times n tables with record counts, where two rows represent genotype levels; n columns 
    represent dependent variable levels. Perform Fisher Exact Tests.
    
    Three tables can be created and three tests can be potentially perform depending on the dataset:
    
    - all records together regardless the gender values, 
    
    - records where gender value is "Male" (if such exists),
    
    - records where gender value is "Female" (if such exists).
    
    Is called from testDataset function when "method" argument is set to "FE" (Fisher Exact Test).
}
\usage{
    buildFisherExactTest(phenList, depVariable, outputMessages=TRUE)
}
\arguments{   
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class}
    \item{depVariable}{a character string defining the dependent variable of interest}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages}   
    
}
\value{
    Return results stored in instance of the \code{\linkS4class{PhenTestResult}} class
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
    result<-testDataset(test,depVariable="Age.In.Weeks",method="FE")
    summaryOutput(result)
}