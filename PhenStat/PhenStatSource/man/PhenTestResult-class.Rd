\name{PhenTestResult-class}
\docType{class}

\alias{PhenTestResult-class}

\title{Class "PhenTestResult"}

\description{A list-based S4 class for storing of model fitting results and other information like dependent variable, method used to build a model, etc.
    Helps to support operations needed for model fitting process and modelling/testing results storage. 
    In the package PhenTestResult object is created by function \code{\link{testDataset}} which in turn calls internal function \code{\link{PhenTestResult}}.}

\section{Explore PhenTestResult object}{
    PhenTestResults object instance contains the following sections:
    
    1. "model.dataset" contains cleaned and checked dataset that was used for the modelling/testing:
    
    \code{phenTestResult$model.dataset}
    
    2. "depVariable" contains dependent variable that was tested: \code{phenTestResult$depVariable}
    
    3. "refGenotype" contains reference genotype which is usually is wildtype: \code{phenTestResult$refGenotype}
    
    4. "method" contains method name that was used for modelling/testing: \code{phenTestResult$method}
    
    5. "equation" is specific for the Mixed Model framework and contains equation that was used for modelling ("withWeight" or "withoutWeight"):
    
    \code{phenTestResult$equation}
    
    6. "pThreshold" contains p-value threshold with default value set to 0.05: 
    
    \code{phenTestResult$pThreshold}
    
    7. "numberSexes" contains the number of sexes in the dataset: \code{phenTestResult$numberSexes}
    
    
    Modelling/testing results are stored in the sections "model.output", "model.formula" and "model.effect". The last two are specific for Mixed Model framework.
    
    The contents of the result sections depend on framework that was used.
    
    
    The results of Mixed Models framework:
    
    1. Batch effect significance: \code{phenTestResult$model.effect.batch}    
    
    2. Variance equality: \code{phenTestResult$model.effect.variance} 
    
    3. Weight effect significance: \code{phenTestResult$model.effect.weight} 
    
    4. Sex effect significance: \code{phenTestResult$model.effect.interaction}     
    
    5. Evidence of sex effect (p-value): \code{phenTestResult$model.output.interaction}
    
    6. Evidence of genotype effect (p-value): 
    \code{phenTestResult$model.output.genotype.nulltest.pVal}
    
    7. Formula for the final genotype model: \code{phenTestResult$model.formula.genotype}
    
    8. Formula for the final null model: \code{phenTestResult$model.formula.null}
    
    9. Model fitting output: \code{phenTestResult$model.output}
    
    10. Model fitting summary: \code{summary(phenTestResult$model.output)$tTable}
    
    
    The results of Fisher Exact Test framework.
    
    Fisher Exact Test framework output contains matrices, statistics about matrices and Fisher Exact Tests results.
    
    1. Count matrices.
    
        1.1 Count matrix for combined dataset: \code{phenTestResult$model.output$count_matrix_all}
    
        1.2 Count matrix for males only dataset: \code{phenTestResult$model.output$count_matrix_male}
    
        1.3 Count matrix for female only dataset: \code{phenTestResult$model.output$count_matrix_female}
    
    2. Percentage matrices.
    
        2.1 Percentage matrix for combined dataset:
     
        \code{phenTestResult$model.output$percentage_matrix_all}
    
        2.2 Percentage matrix for males only dataset: 
    
        \code{phenTestResult$model.output$percentage_matrix_male}
    
        2.3 Percentage matrix for females only dataset: 
    
        \code{phenTestResult$model.output$percentage_matrix_female}

    3. Statistics.
    
        3.1 Statistics for matrix for combined dataset: \code{phenTestResult$model.output$stat_all}
    
        3.2 Statistics for matrix for males only dataset: \code{phenTestResult$model.output$stat_male}
    
        3.3 Statistics for matrix for females only dataset: \code{phenTestResult$model.output$stat_female}

    4. Effect sizes.
    
        4.1 Effect size for combined dataset: \code{phenTestResult$model.output$ES_all}
    
        4.2 Effect size for males only dataset: \code{phenTestResult$model.output$ES_male}
    
        4.3 Effect size for females only dataset: \code{phenTestResult$model.output$ES_female}
 
    5. Fisher Exact Test results. 
    
        5.1 Results for combined dataset: \code{phenTestResult$model.output$all}
    
        5.2 Results for males only dataset: \code{phenTestResult$model.output$male}
    
        5.3 Results for females only dataset: \code{phenTestResult$model.output$female}
    
        5.4 p-value for combined dataset: \code{phenTestResult$model.output$all$p.value} 
    
        5.5 p-value for males only dataset: \code{phenTestResult$model.output$male$p.value}
    
        5.6 p-value for females only dataset: \code{phenTestResult$model.output$female$p.value}
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\seealso{\code{\link{FisherExactTest}}, \code{\link{startModel}}, \code{\link{finalModel}} and \code{\link{PhenTestResult}}}  
 
\examples{
    # Mixed Model framework
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,
            equation="withoutWeight",
            depVariable="Bone.Area")
    result$depVariable
    result$method
    result$equation
    result$numberSexes   
    # Batch effect is significant
    result$model.effect.batch    
    # Variance homogeneosity 
    result$model.effect.variance
    # Weight effect is significant
    result$model.effect.weight     
    # Sex effect is significant 
    result$model.effect.interaction    
    # Sex effect p-value - the result of the test 
    result$model.output.interaction
    # Genotype effect p-value    
    result$model.output.genotype.nulltest.pVal
    # Final model formula with genotype
    result$model.formula.genotype 
    # Final model formula without genotype    
    result$model.formula.null 
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
    result$depVariable
    result$method
    result$numberSexes   
    # Count matrix for all data
    result$model.output$count_matrix_all 
    # Count matrix for males only records
    result$model.output$count_matrix_male
    # Count matrix for females only records
    # result$model.output$count_matrix_female 
    # Percentage matrix for all data
    result$model.output$percentage_matrix_all 
    # Percentage matrix for males only records
    # result$model.output$percentage_matrix_male 
    # Percentage matrix for females only records
    result$model.output$percentage_matrix_female
    # Matrix statistics for all data
    result$model.output$stat_all
    # Matrix statistics for males only records
    result$model.output$stat_male 
    # Matrix statistics for females only records
    result$model.output$stat_female
    # Effect size for all data
    result$model.output$ES_all
    # Effect size for males only records
    result$model.output$ES_male
    # Effect size for females only records
    result$model.output$ES_female 
    # Fisher Exact Test results for all data
    result$model.output$all
    # p-value for all data
    result$model.output$all$p.value
}


\keyword{classes}