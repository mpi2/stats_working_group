\name{performTransformation}
\alias{performTransformation}
\title{Method "performTransformation"}
\description{
    Transforms the given vector of values according to the lambda value:
    log transformation if the lambda is 0, power transformation otherwise
}
\usage{
    performTransformation(values, lambda, scaleShift)
}
\arguments{
    \item{values}{vector of values to reverse}
    \item{lambda}{lambda value of Box-Cox transformation}
    \item{scaleShift}{sclae shift value of Box-Cox transformation}
}
\value{
    Returns values transformed according to Box-Cox rules.
}
\author{Natalja Kurbatova}
\examples{
    performTransformation(c(0.1,38,0.5),-0.7,0)
}    
\keyword{internal}