\name{scatterplotGenotypeWeight}
\alias{scatterplotGenotypeWeight}
\title{Creates a scatter plot body weight versus dependant variable}
\description{
    Creates a scatter plot body weight versus dependant variable
}
\usage{
    scatterplotGenotypeWeight(phenList,depVariable)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class}
    \item{depVariable}{a character string defining the depening variable of interest}
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova}
\seealso{\code{\linkS4class{PhenList}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    scatterplotGenotypeWeight(test,depVariable="Bone.Area")
}