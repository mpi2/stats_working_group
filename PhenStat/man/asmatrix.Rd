\name{as.matrix}
\alias{as.matrix.PhenList}
\title{Method from limma package "as.matrix"}
\description{
    Turn a phenotypic data object into a numeric matrix by extracting the dataset or model output values.
}
\usage{
    \method{as.matrix}{PhenList}(x,\dots)
}
\arguments{
    \item{x}{an object of class \code{PhenList}}
    \item{\dots}{additional arguments, not used for these methods}
}
\details{
    This method extracts the matrix of phenotypic data.
    
    This involves loss of information, so the original data object is not recoverable.
}
\value{
    A numeric matrix.
}
\author{Gordon Smyth}
\seealso{
    \code{\link{as.matrix}} in the base package.
}

\keyword{internal}