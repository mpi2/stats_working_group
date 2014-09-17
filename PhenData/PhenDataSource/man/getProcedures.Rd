\name{getProcedures}
\alias{getProcedures}
\title{Method "getProcedures"}
\description{
    Returns the list of procedures (IMPReSS IDs) that are run for a specified 
    phenotyping center and pipeline.
}
\usage{
    getProcedures(PhenCenterName=NULL,PipelineID=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
}
\value{Returns the list of IMPC procedures (IDS of procedures) that are run by 
specified phenotyping center and pipeline.}

\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    listProcedures  <- getProcedures("WTSI","MGP_001")
    for (procedureIndex in 1:length(listProcedures)) {
        print(paste(listProcedures[procedureIndex],"-",
                        getName("procedure_stable_id","procedure_name",
                        listProcedures[procedureIndex])))
    }
}
