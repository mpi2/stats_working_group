\name{getIMPCDataset}
\alias{getIMPCDataset}
\title{Method "getIMPCDataset"}
\description{
    Returns the IMPC dataset for a knockout line for a specified combination of parameter, procedure, pipeline, 
    phenotyping centre, allele and strain (genetic background) 
    that is used within the data analysis pipeline implemented on the IMPC database.
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
    Returns the IMPC dataset as a data frame.
    - metadata_group    string    a string indicating a group of experimental and control mice that have the same metadata (see also biological_sample_group)
    
    
    - metadata    list of strings    list showing all relevant metadata in effect when the data was collected 
}

\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    df <- getIMPCDataset("WTSI","MGP_001","IMPC_CBC_001","IMPC_CBC_008_001","MGI:4433526")
}
