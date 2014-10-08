\name{printParameters}
\alias{printParameters}
\title{Method "printParameters"}
\description{
    Prints out the list of parameters (IMPReSS IDs and names) that are measured within specified procedure 
    for a pipeline run by phenotyping center.
}
\usage{
    printParameters(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL,n=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
    \item{n}{print only first n parameters; optional argument}
}

\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    printParameters("WTSI","MGP_001","IMPC_CBC_001")
}
