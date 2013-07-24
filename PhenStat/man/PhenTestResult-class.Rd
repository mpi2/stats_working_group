\name{PhenTestResult-class}
\docType{class}

\alias{PhenTestResult-class}
\alias{buildFinalModel,PhenList,PhenTestResult,character,character,numerical,logical vector-method}
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

\description{A list-based S4 class for storing of model fitting results and support of basic operations: tests and model fitting quality}

\section{Create PhenList object}{
    The appropriate way to create PhenTestResult object is via \code{\link{testDataset}} function.
}
\section{Dataset results graphs}{
    The appropriate way to check the dataset properties is vie \code{\link{boxplotGenderGenotype}}, 
    \code{\link{boxplotGenderGenotypeBatch}} and \code{\link{scatterplotGenotypeWeight}} functions.
}
\section{Methods}{
    \describe{
        \item{buildFinalModel}{\code{signature(object1 = "PhenList", object2 = "PhenTestResult", equation = "character",depVariable = "character",
                        pThreshold = "numeric",keepList = "logical vector")}: Build final mixed model}
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
    Placeholder for references
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    result <- testDataset(test,equation="withoutWeight",depVariable="Bone.Area")
    result <- buildFinalModel(test,result,equation="withoutWeight")
    summaryOutput(result)
}


\keyword{classes}