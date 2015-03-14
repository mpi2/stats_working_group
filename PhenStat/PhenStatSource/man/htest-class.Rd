\name{htest-class}
\docType{class}
\alias{htest-class}

\title{Class "htest"}

\description{An old S3 class from package stats rewritten into S4 style in order to add functionality needed for PhenStat.}

\section{Explore htest object}{
    A list with class "htest" containing the following components:
    
    1. p.value: the p-value of the test.
    
    2. conf.int: a confidence interval for the odds ratio. Only present in the 2 by 2 case and if argument conf.int = TRUE.
    
    3. estimate: an estimate of the odds ratio. Note that the conditional Maximum Likelihood Estimate (MLE) rather than the unconditional MLE (the sample odds ratio) is used. Only present in the 2 by 2 case.
    
    4. null.value: the odds ratio under the null, or. Only present in the 2 by 2 case.
    
    5. alternative: a character string describing the alternative hypothesis.
    
    6. method: the character string "Fisher's Exact Test for Count Data".
    
    7. data.name: a character string giving the names of the data.
}

\author{Natalja Kurbatova}

\seealso{\code{\link{htestPhenStat-class}}}

\keyword{classes}
\keyword{internal}