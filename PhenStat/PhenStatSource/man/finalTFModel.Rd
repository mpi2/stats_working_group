\name{finalTFModel}
\alias{finalTFModel}
\title{Method "finalTFModel"}
\description{
    This is an internal function run within TF framework. It completes the final stage of the MM framework, which builds the final model and estimates effects.  
    As an internal function, it doesn't include extensive error testing of inputs. Please use cautiously if calling directly.  
    
    Works with \code{\linkS4class{PhenTestResult}} object created by \code{\link{startModel}} function. 
    The creation of MM final model is based on the significance of different fixed effects, depVariable and equation values 
    stored in \code{\linkS4class{PhenTestResult}} object.
}
\usage{
    finalTFModel(phenTestResult, outputMessages = TRUE)
}
\arguments{   
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class that comes from the function \code{\link{testDataset}}; mandatory argument}
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
\seealso{\code{\linkS4class{PhenTestResult}} and \code{\link{testDataset}}}
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
    result  <- testDataset(test_TF,
                           depVariable="Cholesterol",
                           callAll=FALSE,
                           method="TF")
    # print out formula that has been created
    result$model.formula.genotype
    # print out batch effect's significance 
    result$model.effect.batch
    result <- finalTFModel(result)
}