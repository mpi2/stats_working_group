\name{PhenTestResult-class}
\docType{class}

\alias{PhenTestResult-class}
\alias{buildFinalModel,PhenList,PhenTestResult,character,character,numerical-method}
\alias{testFinalModel,PhenList,PhenTestResult-method}
\alias{qqplotGenotype,PhenList,PhenTestResult-method}
\alias{plotResidualPredicted,PhenList,PhenTestResult-method}
\alias{qqplotRandomEffects,PhenList,PhenTestResult,logical-method}
\alias{boxplotResidualBatch,PhenList,PhenTestResult-method}
\alias{qqplotRotatedResiduals,PhenList,PhenTestResult,logical-method}
\alias{parserOutputSummary,PhenTestResult-method}
\alias{simpleOutput,PhenTestResult-method}
\alias{vectorOutput,PhenTestResult-method}
\alias{outputLength,PhenTestResult-method}

\title{Class "PhenTestResult"}

\description{A list-based S4 class for storing of model fitting results and other information like dependent variable, method used to build a model, etc.
    Help to support operations needed for model fitting process. }

\section{Create PhenTestResult object}{
    The appropriate way to create PhenTestResult object is via \code{\link{testDataset}} function.
}
\section{Explore PhenTestResultObject}{
    Modelling/testing results are stored in the object section model.output. The contents depend on framework that was used.
    The results of MM framework:
    To get was batch effect significant or not: \code{phenTestResult$model.effect.batch}    
    To get was variance equal or not: \code{phenTestResult$model.effect.variance} 
    To get was weight effect significant or not: \code{phenTestResult$model.effect.weight} 
    To get was gender effect significant or not: \code{phenTestResult$model.effect.interaction}     
    To get the evidence of gender effect: \code{phenTestResult$model.output.interaction}
    To get genotype effect: \code{phenTestResult$model.output.genotype.nulltest.pVal}
    
    To get the final model: \code{phenTestResult$model.formula.genotype}
    To get the final null model: \code{phenTestResult$model.formula.null}
    
    To get the model fitting output: \code{phenTestResult$model.output}
    To get the model fitting summary: \code{summary(phenTestResult$model.output)$tTable}
    
    The results of FET framework.
    FET framework output contains matrices, statistics about matrices and Fisher Exact Tests results.
    Count matrices.
    To get count matrix for all data regarding dependent variable: \code{phenTestResult$model.output$count_matrix_all}
    To get count matrix for male only data regarding dependent variable: \code{phenTestResult$model.output$count_matrix_male}
    To get count matrix for female only data regarding dependent variable: \code{phenTestResult$model.output$count_matrix_female}
    
    Percentage matrices.
    To get percentage matrix for all data regarding dependent variable: \code{phenTestResult$model.output$percentage_matrix_all}
    To get percentage matrix for male only data regarding dependent variable: \code{phenTestResult$model.output$percentage_matrix_male}
    To get percentage matrix for female only data regarding dependent variable: \code{phenTestResult$model.output$percentage_matrix_female}

    Statistics.
    To get statistics for matrix for all data regarding dependent variable: \code{phenTestResult$model.output$stat_all}
    To get statistics for matrix for male only data regarding dependent variable: \code{phenTestResult$model.output$stat_male}
    To get statistics for matrix for female only data regarding dependent variable: \code{phenTestResult$model.output$stat_female}

    Effect sizes.
    To get effect size for all data regarding dependent variable: \code{phenTestResult$model.output$ES_all}
    To get effect size for male only data regarding dependent variable: \code{phenTestResult$model.output$ES_male}
    To get effect size for female only data regarding dependent variable: \code{phenTestResult$model.output$ES_female}
 
    Fisher Exact Test results. 
    To get test results for all data regarding dependent variable: \code{phenTestResult$model.output$all}
    To get test results for male only data regarding dependent variable: \code{phenTestResult$model.output$male}
    To get test results for female only data regarding dependent variable: \code{phenTestResult$model.output$female}
    To get p-value only for all data regarding dependent variable: \code{phenTestResult$model.output$all$p.value} 
    To get p-value only for male only records regarding dependent variable: \code{phenTestResult$model.output$male$p.value}
    To get p-value only for female only records regarding dependent variable: \code{phenTestResult$model.output$female$p.value}
}
\section{Dataset results graphs}{
    The appropriate way to check the dataset properties is vie \code{\link{boxplotGenderGenotype}}, 
    \code{\link{boxplotGenderGenotypeBatch}} and \code{\link{scatterplotGenotypeWeight}} functions.
}
\section{Methods}{
    \describe{
        \item{buildFinalModel}{\code{signature(object1 = "PhenList", object2 = "PhenTestResult", equation = "character",depVariable = "character",
                        pThreshold = "numeric")}: Build final mixed model}
        \item{testFinalModel}{\code{signature(object1 = "PhenList", object2 = "PhenTestResult")}: Diagnostic test for mixed model quality of fit}
        \item{qqplotGenotype}{\code{signature(object1 = "PhenList", object2 = "PhenTestResult")}: Creates a Q-Q plot of residuals for each genotype}
        \item{plotResidualPredicted}{\code{signature(object1 = "PhenList", object2 = "PhenTestResult")}: Creates predicted versus residual values plots split by genotype}
        \item{qqplotRandomEffects}{\code{signature(object1 = "PhenList", object2 = "PhenTestResult",keep_batch = "logical")}: Creates a Q-Q plot of blups}
        \item{boxplotResidualBatch}{\code{signature(object1 = "PhenList", object2 = "PhenTestResult")}: Creates a box plot with residue versus batch split by genotype}
        \item{qqplotRotatedResiduals}{\code{signature(object1 = "PhenList", object2 = "PhenTestResult",keep_batch = "logical")}: Creates a Q-Q plot of rotated residuals}
        \item{parserOutputSummary}{\code{signature(object = "PhenTestResult")}: Parsers model output summary and returns it in readable named vector format}
        \item{simpleOutput}{\code{signature(object = "PhenTestResult")}: Returns model fitting and testing results in a user friendly format}
        \item{vectorOutput}{\code{signature(object = "PhenTestResult")}: Returns model fitting and testing results in a vector format}
        \item{outputLength}{\code{signature(object = "PhenTestResult")}: Calculates table length based on used model}        
    }
}
\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,equation="withoutWeight",depVariable="Bone.Area")
    result$model.effect.batch # Batch effect is significant
    result$model.effect.variance # Variance homogeneosity
    result$model.effect.weight # Weight effect is significant
    result$model.effect.interaction # Gender effect is significant 
    result$model.output.interaction # Gender effect p-value - the result of the test
    result$model.output.genotype.nulltest.pVal # Genotype effect p-value    
    result$model.formula.genotype # Final model formula with genotype
    result$model.formula.null # Final model formula without genotype    
    result$model.output # Final model fitting output
    summary(result$model.output)$tTable # Final model fitting summary
    
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/categorical_example1_aff3_Xray.csv"),
            testGenotype="Aff3/Aff3")
    result <- testDataset(test,depVariable="Thoracic.Processes",method="FE")  
    
    result$model.output$count_matrix_all # Count matrix for all data
    result$model.output$count_matrix_male # Count matrix for male only records
    result$model.output$count_matrix_female # Count matrix for female only records
    result$model.output$percentage_matrix_all # Percentage matrix for all data
    result$model.output$percentage_matrix_male # Percentage matrix for male only records
    result$model.output$percentage_matrix_female # Percentage matrix for female only records
    result$model.output$stat_all # Matrix statistics for all data
    result$model.output$stat_male # Matrix statistics for male only records
    result$model.output$stat_female # Matrix statistics for female only records
    result$model.output$ES_all # Effect size for all data
    result$model.output$ES_male # Effect size for male only records
    result$model.output$ES_female # Effect size for female only records
    result$model.output$all #Fisher Exact Test results for all data
    result$model.output$all$p.value #p-value for all data
}


\keyword{classes}