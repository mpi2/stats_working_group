\name{FisherExactTest}
\alias{FisherExactTest}
\title{Method "FisherExactTest"}
\description{
    The main function of the Fisher Exact Test framework. Creates n times 2 matrices with record counts, where n rows 
    represent dependent variable levels and two columns represent genotype levels. Performs Fisher Exact Tests on calculated count matrices.
    
    Three matrices can be created and three tests can be potentially perform depending on the dataset:
    
    - all records together regardless the gender values - combined dataset, 
    
    - records where gender value is "Male" (if such exists) - males only dataset,
    
    - records where gender value is "Female" (if such exists) - females only dataset.
    
    Together with count matrices creates percentage matrices, calculates effect sizes and statistics for count matrices. Performs Fisher Exact Tests. 
    
    The results (matrices, statistics and Fisher Exact Test outputs) are stored in \code{\linkS4class{PhenTestResult}} object.
     
    The funciton is called from "testDataset" function when "method" argument is set to "FE" meaning "Fisher Exact Test".
}
\usage{
    FisherExactTest(phenList, depVariable, outputMessages=TRUE)
}
\arguments{   
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; mandatory argument}
    \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
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
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    file <- system.file("extdata", "test_categorical.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Aff3/Aff3")
    # "FisherExactTest" function is called from "testDataset" function.
    result <- testDataset(test,depVariable="Thoracic.Processes",method="FE")  
    # Fisher Exact Test results can be printed out using function "summaryOutput"
    summaryOutput(result)
}