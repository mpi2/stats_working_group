\name{PhenTestResult-class}
\docType{class}
\alias{PhenTestResult-class}
\alias{getVariable,PhenTestResult-method}
\alias{refGenotype,PhenTestResult-method}
\alias{testGenotype,PhenTestResult-method}
\alias{parameters,PhenTestResult-method}
\alias{method,PhenTestResult-method}
\alias{methodText,PhenTestResult-method}
\alias{analysisResults,PhenTestResult-method}
\alias{analysedDataset,PhenTestResult-method}
\alias{transformationText,PhenTestResult-method}
\alias{show,PhenTestResult-method}
\alias{getCountMatrices,PhenTestResult-method}
\alias{getGenotypeEffect,PhenTestResult-method}
\title{Class "PhenTestResult"}
\description{A list-based S4 class for storing of model fitting results and other information like dependent variable, method used to build a model, etc.
    Helps to support operations needed for model fitting process and modelling/testing results storage. 
    In the package PhenTestResult object is created by function \code{\link{testDataset}}.}
\section{Explore PhenTestResult object}{
    PhenTestResults object instance contains the following slots:
    
    1. "analysedDataset" contains subdataset that was used for the modelling/testing:
    
    \code{analysedDataset(phenTestResult)}
    
    2. "depVariable" contains dependent variable that was tested: \code{getVariable(phenTestResult)}
    
    3. "refGenotype" contains reference genotype which is usually is wildtype: \code{refGenotype(phenTestResult)}
    
    4. "testGenotype" contains test genotype: \code{testGenotype(phenTestResult)} 
    
    5. "method" contains method name that was used for modelling/testing: \code{method(phenTestResult)}
    
    6. "transformationRequired" contains TRUE/FALSE value indicating transofrmation of dependent variable:
    
    \code{transformationText(phenTestResult)} or \code{phenTestResult@transformationRequired}
    
    7. "lambdaValue" contains lambda value for the Box-Cox transformation: 
    
    \code{transformationText(phenTestResult)} or \code{phenTestResult$lambdaValue}
    
    8. "scaleShift" contains the value of scale shift for Box-Cox transformation: 

    \code{transformationText(phenTestResult)} or \code{phenTestResult@scaleShift}
 
    9. "transformationCode" contains the code of transformation. Possible values:
    
    0 - transformation is not applicable (methods "FE", "RR", "LR") or switched off by user
    1 - transformation is not needed (1 is within the 95% confidence interval for lambda)
    2 - log transformation
    3 - power transformation
    4 - transformation is not performed since optimal lambda value is not found  (-5 > lambda > 5)
    
    \code{transformationText(phenTestResult)} or \code{phenTestResult@transformationCode}

    10. "parameters" contains the parameters used during method application, e.g. pThreshold for MM and TF: 
    
    \code{parameters(phenTestResult)}
    
    
    11. Modelling/testing results are stored in the sections "analysisResults": 
    \code{analysisResults(phenTestResult)}
    
    The contents of the analysisResults slot depend on framework that was used.
    
    
    The results of MM and TF frameworks:
    1. Equation used during linear modeling: \code{analysisResults(phenTestResult)$equation}     
    
    2. Batch effect significance: \code{analysisResults(phenTestResult)$model.effect.batch}    
    
    3. Variance equality: \code{analysisResults(phenTestResult)$model.effect.variance} 
    
    4. Weight effect significance: \code{analysisResults(phenTestResult)$model.effect.weight} 
    
    5. Sex effect significance: \code{analysisResults(phenTestResult)$model.effect.interaction}     
    
    6. Evidence of sex effect (p-value): \code{analysisResults(phenTestResult)$model.output.interaction}
    
    7. Evidence of genotype effect (p-value): 
    \code{analysisResults(phenTestResult)$model.output.genotype.nulltest.pVal}
    
    8. Formula for the final genotype model: \code{analysisResults(phenTestResult)$model.formula.genotype}
    
    9. Formula for the final null model: \code{analysisResults(phenTestResult)$model.formula.null}
    
    10. Model fitting output: \code{analysisResults(phenTestResult)$model.output}
    
    11. Model fitting summary: \code{summary(analysisResults(phenTestResult)$model.output)$tTable}
    
    
    The results of FE and RR frameworks are stored as a list of htestPhenStat S4 objects: 
    \code{analysisResults(phenTestResult)[[1]]}.
    
    Each one htestPhenStat object contains:
    
    1. Output of Fisher Exact Test:  \code{pvalue(analysisResults(phenTestResult)[[1]])}.
            
    2. Effect size: \code{ES(analysisResults(phenTestResult)[[1]])}. 
    
    3. Name of the table analysed (all, males, females): \code{analysedSubset(analysisResults(phenTestResult)[[1]])}.
    
    4. Comparison, used for RR only (High vs Normal/Low, Low vs High/Normal): 
            \code{comparison(analysisResults(phenTestResult)[[1]])}.
            
    5. Matrix of counts: \code{matrixCount(analysisResults(phenTestResult)[[1]])}.

}
\section{Slots}{
    \describe{
        \item{\code{analysedDataset}:}{Object of class \code{"data.frame"} ~~ analysed dataset}
        \item{\code{transformationRequired}:}{Object of class \code{"logical"} ~~ flag indicating was or not 
            transformation performed}
        \item{\code{lambdaValue}:}{Object of class \code{"numeric"} ~~ Box-Cox transform lambda value}
        \item{\code{scaleShift}:}{Object of class \code{"numeric"} ~~ Box-Cox transform scale shift}
        \item{\code{transformationCode}:}{Object of class \code{"numeric"} ~~ code explaining the transformation output}
        \item{\code{depVariable}:}{Object of class \code{"character"} ~~ dependent variable}
        \item{\code{refGenotype}:}{Object of class \code{"character"} ~~ reference genotype}
        \item{\code{testGenotype}:}{Object of class \code{"character"} ~~ test genotype}
        \item{\code{method}:}{Object of class \code{"character"} ~~ analysis method}
        \item{\code{parameters}:}{Object of class \code{"matrix"} ~~ parameters of method}
        \item{\code{analysisResults}:}{Object of class \code{"list"} ~~ results of the analysis}
    }
}
\section{Methods}{
    \describe{
        \item{getVariable}{(accessor): Returns dependent variable}
        \item{refGenotype}{(accessor): Returns reference genotype }
        \item{testGenotype}{(accessor): Returns test genotype }
        \item{method}{(accessor): Returns method used }
        \item{methodText}{(accessor): Returns full name of the method used }
        \item{parameters}{(accessor): Returns parameteres used during the method application}
        \item{analysisResults}{(accessor): Returns analysis results depending on method }
        \item{analysedDataset}{(accessor): Returns analysed dataset }
        \item{transformationText}{(accessor): Returns the info about transformation }
        \item{batchIn}{Returns TRUE if the batch is in the dataset, FALSE otherwise}
        \item{weightIn}{Returns TRUE if the weight is in the dataset, FALSE otherwise}
        \item{noSexes}{Returns number of sexes in the dataset}
        \item{show}{Prints out the PhenTestResult object}
        \item{getCountMatrices}{Returns count matrices if they are present (for methods "FE" and "RR"), NULL otherwise.}
    }
}
\author{Natalja Kurbatova}
\examples{
    # Mixed Model framework
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,
            equation="withoutWeight",
            depVariable="Bone.Area")
    getVariable(result)
    method(result)
    # Batch effect is significant
    analysisResults(result)$model.effect.batch    
    # Variance homogeneosity 
    analysisResults(result)$model.effect.variance
    # Weight effect is significant
    analysisResults(result)$model.effect.weight     
    # Sex effect is significant 
    analysisResults(result)$model.effect.interaction    
    # Sex effect p-value - the result of the test 
    analysisResults(result)$model.output.interaction
    # Genotype effect p-value    
    analysisResults(result)$model.output.genotype.nulltest.pVal
    # Final model formula with genotype
    analysisResults(result)$model.formula.genotype 
    # Final model formula without genotype    
    analysisResults(result)$model.formula.null 
    # Final model fitting output
    # result$model.output 
    # Final model fitting summary
    # summary(result$model.output)$tTable 
    
    # Fisher Exact Test framework
    file <- system.file("extdata", "test_categorical.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Aff3/Aff3")
    result <- testDataset(test,
            depVariable="Thoracic.Processes",
            method="FE") 
    getVariable(result)
    method(result)       
    for (i in seq_along(analysisResults(result))) {
        val <- analysisResults(result)[[i]]
        val
    }
}
\keyword{classes}
\keyword{methods}