\name{printHead}
\alias{printHead}
\title{Method from limma package "printHead"}
\description{
    Print leading 5 elements or rows of atomic object
}
\usage{
    printHead(x)
}
\arguments{
    \item{x}{an atomic object}
}
\details{
    Print leading 5 elements or rows of atomic object
    Gordon Smyth
    May 2003.  Last modified 14 April 2009.
}
\references{
    limma package of Bioconductor \url{http://bioconductor.org/packages/release/bioc/html/limma.html}
}
\author{Gordon Smyth}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"),
            testGenotype="Sparc/Sparc")
    printHead(test)
    test
}

\keyword{internal}