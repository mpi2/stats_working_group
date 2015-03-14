\name{formulaDeterminingLambda}
\alias{formulaDeterminingLambda}
\title{Method "formulaDeterminingLambda"}
\description{
    Returns the starting formula that we wish to interogate 
    the data with in estimating the lambda for a Box-Cox transformation
}
\usage{
    formulaDeterminingLambda(noSexes, depVariable, multipleBatches, equation)
}
\arguments{
    \item{noSexes}{number of sexes in the data; mandatory argument}       
    \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
    \item{multipleBatches}{TRUE if there are multiple batch values in the data, FALSE otherwise}      
    \item{equation}{a character string defining the equation to use. 
        Possible values "withWeight" (default),"withoutWeight"}
}
\value{
    Returns a formula
}
\author{Natasha Karp}   
\keyword{internal}