\name{getGenes}
\alias{getGenes}
\title{Method "getGenes"}
\description{
    Returns the list of IMPC genes (IDS of genes) that are processed measuring specified parameter 
    within the procedure in the pipeline run by phenotyping center.
}
\usage{
    getGenes(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL,ParameterID=NULL,StrainID=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
    \item{ParameterID}{IMPC parameter ID; mandatory argument}
    \item{StrainID}{IMPC strain ID; optional argument}
}
\value{
    Returns the list of IMPC genes (IDS of genes) that are processed measuring specified parameter 
    within the procedure in the pipeline run by phenotyping center.
}

\author{Natalja Kurbatova, Jeremy Mason}
\example{
    listGenes  <- getGenes("WTSI","MGP_001","IMPC_CBC_001","IMPC_CBC_003_001")
    for (geneIndex in 1:length(listGenes)) {
        print(paste(listGenes[geneIndex],"-",
                        getName("gene_accession_id","gene_symbol",listGenes[geneIndex])))
    }
}
