\name{getPipelines}
\alias{getPipelines}
\title{Method "getPipelines"}
\description{
    Returns the list of IMPC pipelines (IDS of pipelines) that are run by specified phenotyping center.
}
\usage{
    getPipelines(PhenCenterName=NULL)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
}
\value{
    Returns the list of IMPC pipelines (IDS of pipelines) that are run by specified phenotyping center.
}

\author{Natalja Kurbatova, Jeremy Mason}
\example{
    listPipelines  <- getPipelines("WTSI")
    for (pipelineIndex in 1:length(listPipelines)) {
        print(paste(listPipelines[pipelineIndex],"-",
                    getName("pipeline_stable_id","pipeline_name",listPipelines[pipelineIndex])))
    }
}
