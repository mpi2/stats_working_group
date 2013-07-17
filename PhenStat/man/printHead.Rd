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
\author{Gordon Smyth}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    printHead(test)
    test
}

