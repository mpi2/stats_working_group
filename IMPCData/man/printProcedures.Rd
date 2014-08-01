\name{printProcedures}
\alias{printProcedures}
\title{Method "printProcedures"}
\description{
    Prints out IMPC procedures (ID and name) that are run by specified phenotyping center and pipeline.
}
\usage{
    printProcedures(PhenCenterName=NULL,PipelineID=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
}

\author{Natalja Kurbatova, Jeremy Mason}
\example{
    printProcedures("WTSI","MGP_001")
}
