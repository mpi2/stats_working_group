\name{checkDataset}
\alias{checkDataset}
\title{Method "checkDataset"}
\description{
    Checks dataset for the minimum required info 
}
\usage{
    checkDataset(dataset, outputMessages=TRUE, refGenotype="+/+", dataset.clean=TRUE) 
    
}
\arguments{
    \item{dataset}{data frame created from file or from another source}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages}
    \item{refGenotype}{reference genotype, assigned default value is "+/+"}
    \item{dataset.clean}{flag: "FALSE" value for no modification of the dataset; "TRUE" value to clean dataset if needed}
   
}
\value{
    Returns an instance of the \code{\linkS4class{PhenList}} class.
}


\references{
    Publication placeholder
    
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\seealso{\code{\linkS4class{PhenTestResult}}}

\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    checkDataset(test$dataset)
}