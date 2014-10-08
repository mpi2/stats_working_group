\name{printStrains}
\alias{printStrains}
\title{Method "printStrains"}
\description{
    Prints out the genetic background (MGI IDs and strain name) from which the 
    knockout mice were derived 
    for a specific combination of pipeline, procedure and paramater for a 
    phenotyping center.
}
\usage{
    printStrains(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL,
    ParameterID=NULL,n=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
    \item{ParameterID}{IMPC parameter ID; mandatory argument}
    \item{n}{print only first n strains; optional argument}
}
\value{
    Returns the character containing IMPC strains (IDs and names of strains) 
}
\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    printStrains("WTSI","MGP_001","IMPC_CBC_001","IMPC_CBC_003_001")
}
