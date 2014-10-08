\name{printProcedures}
\alias{printProcedures}
\title{Method "printProcedures"}
\description{
    Prints out the list of procedures (IMPReSS IDs and names) that are run 
    for a specified phenotyping center and pipeline.
}
\usage{
    printProcedures(PhenCenterName=NULL,PipelineID=NULL,n=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{n}{print only first n procedures; optional argument}
}

\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    printProcedures("WTSI","MGP_001")
}
