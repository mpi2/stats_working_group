\name{printAlleles}
\alias{printAlleles}
\title{Method "printAlleles"}
\description{
    Prints out IMPC alleles (ID and name) that are processed measuring specified parameter 
    within the procedure in the pipeline run by phenotyping center.
}
\usage{
    printAlleles(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL,ParameterID=NULL,StrainID=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
    \item{ParameterID}{IMPC parameter ID; mandatory argument}
    \item{StrainID}{IMPC strain ID; optional argument}
}
\author{Natalja Kurbatova, Jeremy Mason}
\example{
    printAlleles("WTSI","MGP_001","IMPC_CBC_001","IMPC_CBC_003_001")
}
