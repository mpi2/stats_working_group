\name{testDataset}
\alias{testDataset}
\title{Method "testDataset"}
\description{      
    Perform the following checks for dependent variable:
    
    1. "depVariable" column should present in the dataset
    
    2. "depVariable" should be numeric for Mixed Model (MM) framework, otherwise recommends Fisher Exact Test (FE)
    
    3. Each one genotype level should have more than one "depVariable" level
    (variability) for MM framework, otherwise recommends FE framework. 
    
    
    Perform other consistency checks. For instance, when "depVariable" column contains cathegorical values 
    but defined method is MM function redefine method to FE.
    
    MM framework: 
    
    For MM framework create start model and modify it after testing of different hypothesis (the model fixed effects).
    
    The model fixed effects are: 
    
    -batch effect (random effect significance), 
    
    -variance effect (TRUE if residual variances for genotype groups are homogeneous and FALSE if they are heterogeneous), 
    
    -interaction effect (genotype by gender interaction significance), 
    
    -gender effect (gender significance), 
    
    -weigth effect (weigth significance).
    
    If user would like to assign other TRUE/FALSE values to the fixed effects of the model 
    then he or she has to define keepList argument which is vector of TRUE/FALSE values.
    
    If user has defined model fixed effects (keepList argument) then function prints out calculated and user defined effects 
    (only when outputMessages argument is set to TRUE), checks user defined effects for consistency 
    (for instance, if there are no "Weight" column in the dataset then weigth effect can't be assigned to TRUE, etc.)
    and modifies start model according to user defined effects. 
    
    As the result PhenTestResult object that contains calculated or user defined model fixed effects and MM start model is created. 

    FE framework:  
    
    For Fisher Exact Test framework there are no special checks except "depVariable" checks mentioned above. 
    "buildFisherExactTest" function is called for count matrix (matrices) creation and Fihser tests performance. 

    As the result PhenTestResult object that contains Fisher tests output is created.   

}    
\usage{
    testDataset(phenList, depVariable, 
            equation="withWeight", outputMessages=TRUE, 
            pThreshold=0.05, method="MM", keepList=NULL)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class}       
    \item{depVariable}{a character string defining the dependent variable of interest}
    \item{equation}{a character string defining the equation to use. Possible values "withWeight" (default),"withoutWeight"}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages}        
    \item{pThreshold}{a numerical value for the p-value threshold to use in tests, default value 0.05}    
    \item{method}{a character string defining the method to use for model building, default value "MM" mixed model}
    \item{keepList}{a logical vector defining the significance of different model effects: keep_batch, keep_equalvar,
        keep_weight, keep_gender, keep_interaction}
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
    result <- testDataset(test,depVariable="Bone.Area")
}