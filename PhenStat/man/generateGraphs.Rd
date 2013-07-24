\name{generateGraphs}
\alias{generateGraphs}
\title{Method "generateGraphs"}
\description{
    Generate different graphs and store them in the defined directory
}
\usage{
    generateGraphs(phenList, phenTestResult, dir, graphingName=NULL, type="Xlib")
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenTestResult}} class}
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class}
    \item{dir}{directory where to store generated graphs}
    \item{graphingName}{a character string defining the name to be used on the graph for the dependent variable}
    \item{type}{character string, one of "Xlib" or "quartz" (some OS X builds) or "cairo"}
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
