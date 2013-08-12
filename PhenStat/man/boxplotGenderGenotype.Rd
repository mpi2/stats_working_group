\name{boxplotGenderGenotype}
\alias{boxplotGenderGenotype}
\title{Method "boxplotGenderGenotype"}
\description{
    Create a box plot split by gender and genotype.
}
\usage{
    boxplotGenderGenotype(phenList,depVariable,graphingName)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; mandatory argument}
    \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
    \item{graphingName}{a character string defining the name to be used on the graph for the dependent variable; mandatory argument}
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenList}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"),
            testGenotype="Sparc/Sparc")
    boxplotGenderGenotype(test,depVariable="Bone.Mineral.Content", 
            graphingName="BMC")
    dataFile <- read.csv("./PhenStat/extdata/MAHN_DEXAdata_1Gender.csv")
    test_1gender <- PhenList(dataset=dataFile,testGenotype="Mysm1/+")
    boxplotGenderGenotype(test_1gender,depVariable="Lean.Mass", 
            graphingName="Lean Mass (g)")
}