\name{printPipelines}
\alias{printPipelines}
\title{Method "printPipelines"}
\description{
    Prints out IMPC pipelines (IDS and names) that are run by specified phenotyping center.
}
\usage{
    printPipelines(PhenCenterName=NULL,excludeLegacyPipelines=TRUE)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{excludeLegacyPipelines}{flag: "FALSE" value to print all pipelines including legacy pipelines; 
        "TRUE" value (default) to print current IMPC pipelines excluding legacy pipelines}
}

\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    printPipelines("WTSI")
}
