\name{PhenTestResult}
\alias{PhenTestResult}
\title{Method "PhenTestResult"}
\description{
    Creates \code{\linkS4class{PhenTestResult}} object based on test results 
    or model building results. 
    Internal function for the package. Not build for users to directly call.
}
\usage{
    PhenTestResult(model.output=NULL, model.dataset=NULL, depVariable=NULL, 
    equation="withWeight", 
            method="MM", model.effect.batch=NULL, 
            model.effect.variance=NULL,model.effect.sex=NULL, 
            model.effect.interaction=NULL, model.output.interaction=NULL, 
            model.effect.weight=NULL, numberSexes=NULL,pThreshold=0.05, 
            model.formula.null=NULL,model.formula.genotype=NULL, 
            model.output.genotype.nulltest.pVal=NULL, 
            model.output.quality=NULL,model.output.summary=NULL,
            model.output.averageRefGenotype = NULL,
            model.output.percentageChanges = NULL) 
}
\arguments{
    \item{model.output}{representing the model fit gls object or lme object}
    \item{model.dataset}{dataset used for modeling}
    \item{depVariable}{depending variable in the model}
    \item{equation}{possible values: "withWeight" to include weight effect 
    into model, "withoutWeight" to exclude weight effect from the model}
    \item{method}{possible values: "MM" to work with Mixed Model framework, 
    "FE" to work with Fisher Exact Test framework}
    \item{model.effect.batch}{Result of the test for batch effect significance}
    \item{model.effect.variance}{Result of the test for residual variance effect}
    \item{model.effect.sex}{Result of the test for sex effect significance}
    \item{model.effect.weight}{Result of the test for weight effect significance}
    \item{model.effect.interaction}{Result of the test for genotype by sex 
    interaction significance}
    \item{model.output.interaction}{Interaction test result: p-value}
    \item{numberSexes}{number of possible sex values in the dataset}
    \item{pThreshold}{p-value threshold used to evaluate the effect 
    significance, default value 0.05}
    \item{model.formula.null}{formula for the model without genotype effect}
    \item{model.formula.genotype}{formula for the model with genotype effect}
    \item{model.output.genotype.nulltest.pVal}{Genotype test results: p-value}
    \item{model.output.quality}{vector that contains diagnostic test output for 
    mixed model quality of fit}
    \item{model.output.summary}{named vector that contains summary of the model 
    output}
    \item{model.output.averageRefGenotype}{named vector that contains mean 
    values calculated for reference genotype records}
    \item{model.output.percentageChanges}{named vector that contains 
    percentage changes values per sex}
}
\value{Returns an instance of the \code{\linkS4class{PhenTestResult}} class.
}
\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and 
    Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} 
    \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide 
    using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\seealso{\code{\linkS4class{PhenTestResult}}}

\examples{
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    # "testDataset" function calls "PhenTestResult" function internally
    result <- testDataset(test,
            equation="withoutWeight",
            depVariable="Lean.Mass")
    class(result)
}