\name{printStrains}
\alias{printStrains}
\title{Method "printStrains"}
\description{
    Prins out IMPC strains (ID and name) that are processed measuring specified parameter 
    within the procedure in the pipeline run by phenotyping center.
}
\usage{
    printStrains(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL,ParameterID=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
    \item{ParameterID}{IMPC parameter ID; mandatory argument}
}

\author{Natalja Kurbatova, Jeremy Mason}
\example{
    printStrains("WTSI","MGP_001","IMPC_CBC_001","IMPC_CBC_003_001")
}
