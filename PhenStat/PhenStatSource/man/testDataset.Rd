\name{testDataset}
\alias{testDataset}
\title{Method "testDataset"}
\description{      
    Performs the following checks for dependent variable:
    
    1. "depVariable" column should present in the dataset
    
    2. "depVariable" should be numeric for Mixed Model (MM) framework, otherwise performs Fisher Exact Test (FE)
    
    3. Each one genotype level should have more than one "depVariable" level
    (variability) for MM framework, otherwise recommends FE framework. 
    
    Runs the stages of selected framework to analyse dependent variable. 
    
    For instance, runs startModel and finalModel for the MM framework if the argument "callAll" is set to TRUE.

    If framework contains only one stage like in Fisher Exact Test case then runs that one stage regardless "callAll" value.

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
    \item{method}{a character string ("MM" or "FE") defining the method to use for model building, default value "MM" mixed model}
    \item{callAll}{flag: "FALSE" value to run only the first stage of the selected framework; "TRUE" value (default) to run all stages of framework one after another}        
    \item{keepList}{a logical vector defining the significance of different model effects: keep_batch, keep_equalvar,
        keep_weight, keep_gender, keep_interaction; default value NULL}
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
    # Mixed Model framework
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,
            depVariable="Lean.Mass")
    # print out formula that has been created
    result$model.formula.genotype
    summaryOutput(result)
    
    # Mixed Model framework with user defined effects
    user_defined_effects <- c(keep_batch=TRUE,
            keep_equalvar=TRUE,
            keep_weight=TRUE,
            keep_gender=TRUE,
            keep_interaction=TRUE)
    result3 <- testDataset(test,
            depVariable="Lean.Mass",
            keepList=user_defined_effects)
    # print out formula that has been created
    result3$model.formula.genotype
    summaryOutput(result3)
    
    # Fisher Exact Test framework
    file <- system.file("extdata", "test_categorical.csv", package="PhenStat")
    test2 <- PhenList(dataset=read.csv(file),
            testGenotype="Aff3/Aff3")
    result2 <- testDataset(test2,
            depVariable="Thoracic.Processes",
            method="FE") 
    summaryOutput(result2) 
}