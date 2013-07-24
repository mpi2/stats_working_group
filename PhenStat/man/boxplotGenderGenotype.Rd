\name{boxplotGenderGenotype}
\alias{boxplotGenderGenotype}
\title{Method "boxplotGenderGenotype"}
\description{
    Creates a box plot split by gender and genotype 
}
\usage{
    boxplotGenderGenotype(phenList,depVariable)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class}
    \item{depVariable}{a character string defining the dependent variable of interest}
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenList}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    boxplotGenderGenotype(test,depVariable="Bone.Area")
}