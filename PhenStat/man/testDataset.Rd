\name{testDataset}
\alias{testDataset}
\title{Method "testDataset"}
\description{      
    Perform the following checks for dependent variable:
    
    1. "depVariable" column should present in the dataset
    
    2. "depVariable" should be numeric for Mixed Model (MM) framework, otherwise performs Fisher Exact Test (FE)
    
    3. Each one genotype level should have more than one "depVariable" level
    (variability) for MM framework, otherwise recommends FE framework. 


}    
\usage{
    testDataset(phenList, depVariable, 
            equation="withWeight", outputMessages=TRUE, 
            pThreshold=0.05, method="MM", callAll=TRUE, keepList=NULL)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; mandatory argument}       
    \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
    \item{equation}{a character string defining the equation to use. Possible values "withWeight" (default),"withoutWeight"}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages; default value TRUE}        
    \item{pThreshold}{a numerical value for the p-value threshold to use in tests, default value 0.05}    
    \item{method}{a character string defining the method to use for model building, default value "MM" mixed model}
    \item{callAll}{flag: "FALSE" value to run only the first method (framework) function; "TRUE" value (default) to run all methods (framework) functions}        
    \item{keepList}{a logical vector defining the significance of different model effects: keep_batch, keep_equalvar,
        keep_weight, keep_gender, keep_interaction; default value NULL}
}
\value{
    Return results stored in instance of the \code{\linkS4class{PhenTestResult}} class
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenList}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,depVariable="Lean.Mass")
}