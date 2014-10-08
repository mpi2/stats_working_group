\name{printPipelines}
\alias{printPipelines}
\title{Method "printPipelines"}
\description{
    Prints out the list of IMPC pipelines (IMPReSS IDs and names) that are run 
    by specified phenotyping center.
}
\usage{
    printPipelines(PhenCenterName=NULL,n=NULL, excludeLegacyPipelines=TRUE)
}
\arguments{
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{n}{print only first n pipelines; optional argument}
    \item{excludeLegacyPipelines}{flag: "FALSE" value to print all pipelines 
    including legacy pipelines; 
        "TRUE" value (default) to print current IMPC pipelines excluding legacy 
        pipelines}
}
\value{
    Returns the character containing IMPC pipelines (IDs and names of pipelines) 
}
\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    printPipelines("WTSI")
}
