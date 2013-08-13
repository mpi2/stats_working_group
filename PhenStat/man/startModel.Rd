\name{startModel}
\alias{startModel}
\title{Method "startModel"}
\description{      
    The first function to be run within MM framework. 
    
    It creates start model and modify it after testing of different hypothesis (the model fixed effects).
    
    The model fixed effects are: 
    
    -batch effect (random effect significance), 
    
    -variance effect (TRUE if residual variances for genotype groups are homogeneous and FALSE if they are heterogeneous), 
    
    -interaction effect (genotype by gender interaction significance), 
    
    -gender effect (gender significance), 
    
    -weight effect (weight significance).
    
    If user would like to assign other TRUE/FALSE values to the fixed effects of the model 
    then he or she has to define keepList argument which is vector of TRUE/FALSE values.
    
    If user has defined model fixed effects (keepList argument) then function prints out calculated and user defined effects 
    (only when outputMessages argument is set to TRUE), checks user defined effects for consistency 
    (for instance, if there are no "Weight" column in the dataset then weight effect can't be assigned to TRUE, etc.)
    and modifies start model according to user defined effects. 
    
    As the result PhenTestResult object that contains calculated or user defined model fixed effects and MM start model is created.     
}    
\usage{
    startModel(phenList, depVariable, 
            equation="withWeight", outputMessages=TRUE, 
            pThreshold=0.05, keepList=NULL)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; mandatory argument}       
    \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
    \item{equation}{a character string defining the equation to use. Possible values "withWeight" (default),"withoutWeight"}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages; default value TRUE}        
    \item{pThreshold}{a numerical value for the p-value threshold to use in tests, default value 0.05}    
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
    result <- testDataset(test,depVariable="Lean.Mass",callAll=FALSE)
}