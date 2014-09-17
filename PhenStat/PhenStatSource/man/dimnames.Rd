\name{dimnames}
\alias{dimnames.PhenList}
\alias{dimnames<-.PhenList}
\title{Method "dimnames"}
\description{
    Retrieve the dimension names of a phenotypic data object.
    Based on dimnames methods of LIMMA package.
}
\usage{
    \method{dimnames}{PhenList}(x)
    \method{dimnames}{PhenList}(x) <- value
}
\arguments{
    \item{x}{an object of class \code{PhenList}}
    \item{value}{a possible value for \code{dimnames(x)}: 
    see \code{\link{dimnames}}}
}
\details{
    The dimension names of a phenotypic data object are the same as those of 
    the most important matrix component of that object.
    
    A consequence is that \code{rownames} and \code{colnames} will work as 
    expected.
}
\value{Either \code{NULL} or a list of length 2.
    If a list, its components are either \code{NULL} or a character vector the 
    length of the appropriate dimension of \code{x}.
}
\author{Natalja Kurbatova}
\seealso{
    \code{\link[base]{dimnames}} in the base package.
}

\keyword{internal}