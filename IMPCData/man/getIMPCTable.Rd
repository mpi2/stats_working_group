\name{getIMPCTable}
\alias{getIMPCTable}
\title{Method "getIMPCTable"}
\description{
    Returns the IMPC dataset.
}
\usage{
    getIMPCTable(fileName="ImpressData.csv",PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, ParameterID=NULL, 
            AlleleID=NULL, StrainID=NULL)
}
\arguments{
    \item{fileName}{name of the file where to save resulting table with IMPC objects; mandatory argument}
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
    getIMPCTable("./IMPCData.csv","WTSI","MGP_001","IMPC_CBC_001")
}
