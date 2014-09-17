\name{getPhenCenters}
\alias{getPhenCenters}
\title{Method "getPhenCenters"}
\description{
    Returns the list of IMPC phenotyping centers. 
}
\usage{
    getPhenCenters(excludeLegacyPipelines=TRUE)
}
\arguments{
    \item{excludeLegacyPipelines}{flag: FALSE value to get all IMPC phenotyping 
    centers that run all kind of pipelines including legacy pipelines; 
    TRUE value (default) to get IMPC phenotyping centers 
        excluding centers that run only legacy pipelines}
}

\value{Returns the list of IMPC phenotyping centers.}

\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    getPhenCenters()
}

