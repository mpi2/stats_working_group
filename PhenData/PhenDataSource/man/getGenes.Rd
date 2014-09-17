\name{getGenes}
\alias{getGenes}
\title{Method "getGenes"}
\description{
    Returns the genes (MGI IDs) that are reported for a specified combination 
    of parameter, procedure, pipeline and phenotyping center.
}
\usage{
    getGenes(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL,
    ParameterID=NULL,StrainID=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
    \item{ParameterID}{IMPC parameter ID; mandatory argument}
    \item{StrainID}{IMPC strain ID; optional argument}
}
\value{Returns the list of IMPC genes (IDS of genes) that are processed 
measuring specified parameter within the procedure in the pipeline run by 
phenotyping center.}

\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    listGenes  <- getGenes("WTSI","MGP_001","IMPC_CBC_001","IMPC_CBC_003_001")
}
