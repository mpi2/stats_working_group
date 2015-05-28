\name{performReverseTransformation}
\alias{performReverseTransformation}
\title{Method "performReverseTransformation"}
\description{
    Reverse back the transformed values according to the lambda value:
    exponential transformation if the lambda is 0, fractional power transformation otherwise 
}
\usage{
    performReverseTransformation(values, lambda, scaleShift)
}
\arguments{
    \item{values}{vector of values to reverse}
    \item{lambda}{lambda value of Box-Cox transformation}
    \item{scaleShift}{sclae shift value of Box-Cox transformation}
}
\value{
    Returns values transformed back to original scale. 
}
\author{Natalja Kurbatova}
\examples{
    performReverseTransformation(c(-5.7312462,1.3166139,-0.8921497),-0.7,0)
}    
\keyword{internal}