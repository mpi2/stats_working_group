\name{startTFModel}
\alias{startTFModel}
\title{Method "startTFModel"}
\description{      
    This is an internal function run within TF framework. It completes the testing stage of which effects are significant.  
    As an internal function, it doesn't include extensive error testing of inputs. Please use cautiously if calling directly.  
    
    It creates start model and modify it after testing of different hypothesis.
    
    The tested fixed effects are: 
    
    -batch effect (TRUE if batch variation is significant, FALSE if not), 
    
    -variance effect (TRUE if residual variances for genotype groups are homogeneous and FALSE if they are heterogeneous), 
    
    -interaction effect (TRUE if genotype by sex interaction is significant), 
    
    -sex effect (TRUE if sex is significant), 
    
    -weight effect (TRUE if weight is significant).
    
}    
\usage{
    startTFModel(phenList, depVariable, 
            equation="withWeight", outputMessages=TRUE, 
            pThreshold=0.05, keepList=NULL)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; mandatory argument}       
    \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
    \item{equation}{a character string defining the equation to use. 
        Possible values "withWeight" (default),"withoutWeight"}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages
        ; default value TRUE}        
    \item{pThreshold}{a numerical value for the p-value threshold used to determine which fixed effects to keep 
        in the model, default value 0.05}    
    \item{keepList}{a logical vector defining the significance of different model effects: keep_batch, keep_equalvar,
        keep_weight, keep_sex, keep_interaction; default value NULL}
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
    file <- system.file("extdata", "test7_TFE.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
                     testGenotype="het",
                     refGenotype = "WT",
                     dataset.colname.sex="sex",
                     dataset.colname.genotype="Genotype",
                     dataset.values.female="f",
                     dataset.values.male= "m",
                     dataset.colname.weight="body.weight",
                     dataset.colname.batch="Date_of_procedure_start")

    test_TF <- TFDataset(test,depVariable="Cholesterol")
    
    # when "testDataset" function's argument "callAll" is set to FALSE 
    # only "startTFModel" function is called - the first step of TFE framework
    result <- testDataset(test_TF,
            depVariable="Cholesterol",
            callAll=FALSE,
            method="TF")
    # print out formula that has been created
    result$model.formula.genotype
    # print out batch effect's significance 
    result$model.effect.batch
}