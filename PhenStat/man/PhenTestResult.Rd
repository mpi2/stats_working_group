\name{PhenTestResult}
\alias{PhenTestResult}
\title{Method "PhenTestResult"}
\description{
    Creates PhenTestResult object, based on test results or model building results. 
    Returns an instance of the PhenTestResult object.
}
\usage{
    PhenTestResult(model.output=NULL, depVariable=NULL, equation="withWeight", model.effect.batch=NULL, 
            model.effect.variance=NULL,model.effect.gender=NULL, model.effect.interaction=NULL, 
            model.output.interaction=NULL, model.effect.weight=NULL, numberGenders=NULL,pThreshold=0.05, 
            model.formula.null=NULL,model.formula.genotype=NULL, model.output.genotype.nulltest.pVal=NULL, 
            model.output.quality=NULL,model.output.summary=NULL) 
}
\arguments{
    \item{model.output}{representing the model fit gls object or lme object}
    \item{depVariable}{depending varibale in the model}
    \item{equation}{possible values: "withWeight" to include weight effect into model, "withoutWeight" to exclude weight effect from the model}
    \item{model.effect.batch}{Result of the test for batch effect significance}
    \item{model.effect.variance}{Result of the test for residual variance effect}
    \item{model.effect.gender}{Result of the test for gender effect significance}
    \item{model.effect.weight}{Result of the test for weight effect significance}
    \item{model.effect.interaction}{Result of the test for genotype by gender interaction significance}
    \item{model.output.interaction}{Interaction test result: p-value}
    \item{numberGenders}{number of possible gender values in the dataset}
    \item{pThreshold}{p-value threshold used to evaluate the effect significance, default value 0.05}
    \item{model.formula.null}{formula for the model withput genotype effect}
    \item{model.formula.genotype}{formula for the model with genotype effect}
    \item{model.output.genotype.nulltest.pVal}{Genotype test results: p-value}
    \item{model.output.quality}{vector that contains diagnostic test output for mixed model quality of fit}
    \item{model.output.summary}{named vector that contains summary of the model output}
}
\value{
    Returns an instance of the \code{\linkS4class{PhenTestResult}} class.
}
\references{
        Chapter 3. "Linear Mixed Models" ...
        PLOS one manuscript ...
    
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\seealso{\code{\linkS4class{PhenList}}}

\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    result <- testDataset(test,equation="withoutWeight",depVariable="Bone.Area")
}