\name{generateGraphs}
\alias{generateGraphs}
\title{Method "generateGraphs"}
\description{
    Function generates graphs for the Mixed Model and Time as Fixed Effect frameworks' 
    results and stores generated graphs in the defined directory.
}
\usage{
    generateGraphs(phenTestResult, dir, 
            graphingName=NULL, type="Xlib")
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class that comes from the function \code{\link{testDataset}}; mandatory argument}
    \item{dir}{directory where to store generated graphs; mandatory argument}
    \item{graphingName}{a character string defining the name to be used on the graph for the dependent variable; defaults to NULL then the depVariable name will be plotted}
    \item{type}{character string, one of "Xlib" or "quartz" (some OS X builds) or "cairo", or "window", or cairo-png"; mandatory argument}
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
