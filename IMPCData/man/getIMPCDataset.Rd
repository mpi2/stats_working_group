\name{getIMPCDataset}
\alias{getIMPCDataset}
\title{Method "getIMPCDataset"}
\description{
    Returns the IMPC dataset.
}
\usage{
    getIMPCDataset(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL,ParameterID=NULL,
            AlleleID=NULL,StrainID=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
    \item{ParameterID}{IMPC parameter ID; mandatory argument}
    \item{AlleleID}{IMPC allele ID; mandatory argument}
    \item{StrainID}{IMPC strain ID; optional argument}
}
\value{
    Returns the IMPC dataset.
}

\author{Natalja Kurbatova, Jeremy Mason}
\example{
    df2 <- getIMPCDataset('WTSI','MGP_001','IMPC_CBC_001','IMPC_CBC_008_001','MGI:4433526')
}
