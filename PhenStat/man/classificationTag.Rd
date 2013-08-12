\name{classificationTag}
\alias{classificationTag}
\title{Method "classificationTag"}
\description{
    Returns a classification tag to assign a sexual dimorphism assessment of the phenotypic change
}
\usage{
    classificationTag(result, interactionMode = TRUE, 
            phenotypeThreshold = 0.01)
}
\arguments{
    \item{result}{instance of the \code{\linkS4class{PhenTestResult}} class; ; mandatory argument}
    \item{interactionMode}{flag: "FALSE" value for the no interaction messages to the user; "TRUE" value to interact}
    \item{phenotypeThreshold}{a numerical value defining the significance threshold to use in assigning the classification tag}
}
\value{
    Return model fitting results in a text format
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,depVariable="Lean.Mass")
    classificationTag(result)
}