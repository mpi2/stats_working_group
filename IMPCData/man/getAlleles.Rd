\name{getAlleles}
\alias{getAlleles}
\title{Method "getAlleles"}
\description{
    Returns the list of IMPC alleles (IDS of alleles) that are processed measuring specified parameter 
    within the procedure in the pipeline run by phenotyping center.
}
\usage{
    getAlleles(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL,ParameterID=NULL,StrainID=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
    \item{ParameterID}{IMPC parameter ID; mandatory argument}
    \item{StrainID}{IMPC strain ID; optional argument}
}
\value{
    Returns the list of IMPC alleles (IDS of alleles) that are processed measuring specified parameter 
    within the procedure in the pipeline run by phenotyping center.
}

\author{Natalja Kurbatova, Jeremy Mason}
\example{
    listAlleles  <- getAlleles("WTSI","MGP_001","IMPC_CBC_001","IMPC_CBC_003_001")
    for (alleleIndex in 1:length(listAlleles)) {
        print(paste(listAlleles[alleleIndex],"-",
                        getName("allele_accession_id","allele_symbol",listAlleles[alleleIndex])))
    }
}
