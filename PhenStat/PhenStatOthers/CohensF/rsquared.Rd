\name{rsquared}
\alias{rsquared}
\title{Method "rsquared"}
\description{
    Returns R squared coefficient of determination of the model fitting. 
    A measure of how well observed outcomes are replicated by the model. 
    Works for "gls" and "lme" modelling techniques.
}
\usage{
    rsquared(i,base_level_model=NULL)
}
\arguments{
    \item{i}{model to evaluate; mandatory argument}
    \item{base_level_model}{model that indicates the base level. In the case of Mixed Models approach it is intercept only model; needed only for gls technique}
}
\value{
    Returns R squared coefficient of determination of the model fitting. 
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/test.csv"),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,
            depVariable="Lean.Mass")
    R2_genotype_model <- rsquared(phenTestResult$model.output)
    R2_null_model <- rsquared(phenTestResult$model.null)
}