\name{getZygosities}
\alias{getZygosities}
\title{Method "getZygosities"}
\description{
    Returns the list of IMPC zygosities (IDS of genes) of genes/alleles that are processed measuring specified parameter 
    within the procedure in the pipeline run by phenotyping center.
}
\usage{
    getZygosities(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL,ParameterID=NULL,
            StrainID=NULL,GeneID=NULL, AlleleID=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
    \item{ParameterID}{IMPC parameter ID; mandatory argument}
    \item{StrainID}{IMPC strain ID; optional argument}
    \item{GeneID}{IMPC gene ID; optional argument}
    \item{AlleleID}{IMPC allele ID; optional argument}
}
\value{
    Returns the list of IMPC genes (IDS of genes) that are processed measuring specified parameter 
    within the procedure in the pipeline run by phenotyping center.
}

\author{Natalja Kurbatova, Jeremy Mason}
\example{
    listZygosities  <- getZygosities("WTSI","MGP_001","IMPC_CBC_001","IMPC_CBC_003_001",GeneID="MGI:1194894")
    for (zygosityIndex in 1:length(listZygosities)) {
        print(listZygosities[zygosityIndex])
    }
}
