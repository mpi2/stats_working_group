\name{LargeDataObject-class}
\docType{class}

\alias{LargeDataObject-class}
\alias{show,LargeDataObject-method}
\alias{printHead,LargeDataObject-method}

\title{Class "LargeDataObject"}

\description{Internal class. A list-based S4 class for storing large data objects.}

\section{Methods}{
    \describe{
        \item{show}{\code{signature(object = "LargeDataObject")}: Prints out the object head using common method show}
        \item{printHead}{\code{signature(object = "LargeDataObject")}: Print leading 5 elements or rows of atomic object}
    }
}
\references{
    limma package of Bioconductor \url{http://bioconductor.org/packages/release/bioc/html/limma.html}
}
\author{Gordon Smyth}

\keyword{classes}

\keyword{internal}