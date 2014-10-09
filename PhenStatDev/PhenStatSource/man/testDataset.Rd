\name{testDataset}
\alias{testDataset}
\title{Method "testDataset"}
\description{  
    Statistical analysis manager function in PhenStat package. 
    
    Firstly, it performs the checks for dependent variable in the dataset. 
    Some checks are generic, some depends on selected statistical 
    framework/method.
    
    Secondly, if checks are passed it runs the stages of selected framework 
    to analyse dependent variable. 
    
    For instance, runs startModel and finalModel for the MM framework if the 
    argument "callAll" is set to TRUE.
    If framework contains only one stage like in Fisher Exact Test case then 
    runs that one stage regardless "callAll" value.

}    
\usage{
    testDataset(phenList=NULL, depVariable=NULL, equation="withWeight", 
    outputMessages=TRUE,pThreshold=0.05, method="MM", callAll=TRUE, 
    keepList=NULL,dataPointsThreshold=4,RR_naturalVariation=95,
    RR_controlPointsThreshold=60)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; 
    mandatory argument}       
    \item{depVariable}{a character string defining the dependent variable of 
    interest; mandatory argument}
    \item{equation}{a character string defining the equation to use. 
        Possible values "withWeight" (default), "withoutWeight"}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; 
    "TRUE" value to show output messages; 
        default value TRUE}        
    \item{pThreshold}{a numerical value for the p-value threshold used to 
    determine which fixed effects 
        to keep in the model, default value 0.05}    
    \item{method}{a character string ("MM", "FE" or "RR") defining the method
    to use for model building; 
        default value "MM" for mixed model}
    \item{callAll}{flag: "FALSE" value to run only the first stage of the 
    selected framework; "TRUE" value (default) 
        to run all stages of framework one after another}        
    \item{keepList}{a logical vector defining the significance of different
    model effects: keep_batch, keep_equalvar,
        keep_weight, keep_sex, keep_interaction; default value NULL}
    \item{dataPointsThreshold}{threshold for the number of data points in
    the MM framework; default value is 4
        ; minimal value is 2}
    \item{RR_naturalVariation}{threshold for the variation ranges in the
    RR framework; default value is 95
        ; minimal value is 60}
    \item{RR_controlPointsThreshold}{threshold for the number of control
    data points in the RR framework
        ; default value is 60; minimal value is 40}
}
\value{Returns results stored in instance of the \code{\linkS4class{PhenTestResult}} class
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and 
    Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} 
    \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide 
    using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
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
            keep_sex=TRUE,
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