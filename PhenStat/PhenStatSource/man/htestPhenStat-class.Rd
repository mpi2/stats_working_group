\name{htestPhenStat-class}
\docType{class}
\alias{htestPhenStat-class}
\alias{getColumnView,htestPhenStat-method}
\alias{getMatrix,htestPhenStat-method}
\alias{getPercentageMatrix,htestPhenStat-method}
\alias{subsetText,htestPhenStat-method}
\alias{comparison,htestPhenStat-method}
\alias{matrixCount,htestPhenStat-method}
\alias{analysedSubset,htestPhenStat-method}
\alias{pvalue,htestPhenStat-method}
\alias{subsetText,htestPhenStat-method}
\alias{show,htestPhenStat-method}
\title{Class "htestPhenStat"}

\description{Extension of htest class for the PhenStat needs}

\section{Explore htestPhenStat object}{
    A list with "htest" components in modelOutput slot and the additional PhenStat specific slots:
    
    1. modelOutput: htest object with Fisher Exact test (fisher.test) output
    
    2. analysedSubset: the subset analysed with possible values "all", "males", "females"
    
    3. comparison: for the RR method with possible values "High vs Normal/Low", "Low vs Normal/High"
    
    4. ES: effect size value
    
    5. matrixCount: matrix of counts used for the test

}
\section{Slots}{
    \describe{
        \item{\code{modelOutput}:}{Object of class \code{"htest"} ~~ }
        \item{\code{analysedSubset}:}{Object of class \code{"character"} ~~ }
        \item{\code{comparison}:}{Object of class \code{"character"} ~~ }
        \item{\code{ES}:}{Object of class \code{"numeric"} ~~ }
        \item{\code{matrixCount}:}{Object of class \code{"matrix"} ~~ }
    }
}
\section{Methods}{
    \describe{
        \item{getColumnView}{\code{signature(obj = "htestPhenStat")}: ... }
        \item{getMatrix}{\code{signature(obj = "htestPhenStat")}: ... }
        \item{getPercentageMatrix}{\code{signature(obj = "htestPhenStat")}: ... }
        \item{subsetText}{\code{signature(object = "htestPhenStat")}: ... }
        \item{comparison}{\code{signature(object = "htestPhenStat")}: ... }
        \item{show}{\code{signature(object = "htestPhenStat")}: ... }
    }
}

\author{Natalja Kurbatova}

\seealso{\code{\link{htest-class}}}

\examples{
    file <- system.file("extdata", "test_categorical.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Aff3/Aff3")
    # "FisherExactTest" function is called from "testDataset" function.
    result <- testDataset(test,depVariable="Thoracic.Processes",method="FE")  
    # Fisher Exact Test results can be printed out using function "summaryOutput"
    val <- analysisResults(result)[[1]]
    print(val)
    class(val)
}

\keyword{classes}
\keyword{methods}
\keyword{internal}