\name{boxplotGenderGenotypeBatch}
\alias{boxplotGenderGenotypeBatch}
\title{Method "boxplotGenderGenotypeBatch"}
\description{
    Creates a box plot split by gender, genotype and batch
}
\usage{
    boxplotGenderGenotypeBatch(phenList,depVariable)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class}
    \item{depVariable}{a character string defining the dependent variable of interest}
    \item{graphingName}{a character string defining the name to be used on the graph for the dependent variable}
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenList}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    boxplotGenderGenotypeBatch(test,depVariable="Bone.Mineral.Content",graphingName="BMC" )
    test_1gender <- PhenList(dataset=read.csv("./PhenStat/extdata/MAHN_DEXAdata_1Gender.csv"))
    boxplotGenderGenotypeBatch(test_1gender,depVariable="Lean.Mass", graphingName="Lean Mass")
}