\name{Cohenf.Marginal}
\alias{Cohenf.Marginal}
\title{Method "Cohenf.Marginal"}
\description{
    Returns Cohen's f squared value.
    A measure of local effect size calculated as variance among group means divided by pooled within group variance. 
    Works for "gls" and "lme" modelling techniques.
    Marginal here means the variance associated with the fixed effects is considered.
}
\usage{
    Cohenf.Marginal(phenTestResult)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class that comes from the function \code{\link{testDataset}}; mandatory argument}
}
\value{
    Returns Cohen's f squared local effect size. 
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
    Cohenf.Marginal(result)
}