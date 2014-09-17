\name{dim}
\alias{dim.PhenList}
\title{Method "dim"}
\description{
    Retrieve the number of rows (measures) and columns (parameters) for a 
    PhenList object.
}
\usage{
    \method{dim}{PhenList}(x)
}
\arguments{
    \item{x}{an object of class \code{PhenList}}
}
\details{
    These methods allow one to extract the size of phenotypic data objects in 
    the same way that one would do for ordinary matrices.
    
    A consequence is that row and column commands \code{nrow(x)}, \code{ncol(x)} 
    and so on also work.
}
\value{Numeric vector of length 2.
    The first element is the number of rows (measures) and the second is the 
    number of columns (parameters).
}
\author{Natalja Kurbatova}
\seealso{
    \code{\link[base]{dim}} in the base package.
}
\examples{
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    dim(test)
}
\keyword{internal}

