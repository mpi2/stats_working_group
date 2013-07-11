\name{boxplotGenderGenotypeBatch}
\alias{boxplotGenderGenotypeBatch}
\title{Creates a box plot split by gender, genotype and batch}
\description{
    Creates a box plot split by gender, genotype and batch
}
\usage{
    boxplotGenderGenotypeBatch(phenList,depVariable)
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
    boxplotGenderGenotypeBatch(test,depVariable="Bone.Area")
}