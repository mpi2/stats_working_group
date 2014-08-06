\name{printParameters}
\alias{printParameters}
\title{Method "printParameters"}
\description{
    Prints out IMPC parameters (ID and name) that are measured within specified procedure in 
    the pipeline run by phenotyping center.
}
\usage{
    printParameters(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
}

\author{Natalja Kurbatova, Jeremy Mason}
\example{
    printParameters("WTSI","MGP_001","IMPC_CBC_001")
}
