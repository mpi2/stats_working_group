\name{getStrains}
\alias{getStrains}
\title{Method "getStrains"}
\description{
    Returns the genetic background (MGI IDs) from which the knockout mice were derived 
    for a specific combination of pipeline, procedure and paramater for a phenotyping center. 
}
\usage{
    getStrains(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL,ParameterID=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
    \item{ParameterID}{IMPC parameter ID; mandatory argument}
}
\value{
    Returns the list of IMPC strains (IDS of strains) that are processed measuring specified parameter 
    within the procedure in the pipeline run by phenotyping center.
}

\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    listStrains  <- getStrains("WTSI","MGP_001","IMPC_CBC_001","IMPC_CBC_003_001")
    for (strainIndex in 1:length(listStrains)) {
        print(paste(listStrains[strainIndex],"-",
                        getName("strain_accession_id","strain_name",listStrains[strainIndex])))
    }
}
