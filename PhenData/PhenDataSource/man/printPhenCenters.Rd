\name{printPhenCenters}
\alias{printPhenCenters}
\title{Method "printPhenCenters"}
\description{
    Prints out the list of IMPC phenotyping centers. 
}
\usage{
    printPhenCenters(excludeLegacyPipelines=TRUE)
}
\arguments{
    \item{excludeLegacyPipelines}{flag: "FALSE" value to print all IMPC phenotyping centers that run all kind of 
        pipelines including legacy pipelines; "TRUE" value (default) to print IMPC phenotyping centers 
        excluding centers that run only legacy pipelines}
}

\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    printPhenCenters()
}
