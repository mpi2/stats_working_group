\name{printGenes}
\alias{printGenes}
\title{Method "printGenes"}
\description{
    Prints out genes (MGI IDs and symbols) that are reported for a specified 
    combination of parameter, 
    procedure, pipeline and phenotyping center.
}
\usage{
    printGenes(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL,
    ParameterID=NULL,StrainID=NULL,n=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
    \item{ParameterID}{IMPC parameter ID; mandatory argument}
    \item{StrainID}{IMPC strain ID; optional argument}
    \item{n}{print only first n genes; optional argument}
}
\value{
    Returns the character containing IMPC genes (IDs and names of genes) 
}
\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    printGenes("WTSI","MGP_001","IMPC_CBC_001","IMPC_CBC_003_001")
}
