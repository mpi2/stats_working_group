\name{getParameters}
\alias{getParameters}
\title{Method "getParameters"}
\description{
    Returns the list of parameters (IMPReSS IDs) that are measured within specified procedure 
    for a pipeline run by phenotyping center.
}
\usage{
    getParameters(PhenCenterName=NULL,PipelineID=NULL,ProcedureID=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
}
\value{
    Returns the list of IMPC parameters (IDS of parameters) that are measured within specified procedure in 
    the pipeline run by phenotyping center.
}

\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    listParameters  <- getParameters("WTSI","MGP_001","IMPC_CBC_001")
    for (parameterIndex in 1:length(listParameters)) {
        print(paste(listParameters[parameterIndex],"-",
                        getName("parameter_stable_id","parameter_name",listParameters[parameterIndex])))
    }
}
