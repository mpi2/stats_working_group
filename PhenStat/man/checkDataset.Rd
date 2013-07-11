\name{checkDataset}
\alias{checkDataset}

\title{Check dataset for the minimum required info}
\description{
    Check dataset for the minimum required info 
}
\usage{
    checkDataset(dataset, interactionMode=TRUE, refGenotype="+/+", dataset.clean=FALSE) 
    
}
\arguments{
    \item{dataset}{data frame created from file or from another source}
    \item{interactionMode}{flag: "FALSE" value for the no interaction messages to the user; "TRUE" value to interact}
    \item{refGenotype}{reference genotype, assigned default value is "+/+"}
    \item{dataset.clean}{flag: "FALSE" value for no modification of the dataset; "TRUE" value to clean dataset if needed}
   
}
\value{
    Returns an instance of the \code{\linkS4class{PhenList}} class.
}


\references{
    Publication placeholder
    
}
\author{Natalja Kurbatova}

\seealso{\code{\linkS4class{PhenTestResult}}}

\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    checkDataset(test$dataset)
}