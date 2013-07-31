\name{PhenList-class}
\docType{class}
\alias{PhenList-class}
\alias{testDataset,PhenList,character,character,numeric,character-method}
\alias{scatterplotGenotypeWeight,PhenList,character-method}
\alias{boxplotGenderGenotypeBatch,PhenList,character-method}
\alias{boxplotGenderGenotype,PhenList,character-method}  

\title{Class "PhenList"}

\description{A list-based S4 class for storing phenotypic data. Help to support data integrity checks and statistics calculation.}

\section{Create PhenList object}{
    The appropriate way to create PhenList object is via \code{\link{PhenList}} function.
}
\section{Dataset graphs}{
    The appropriate way to check the dataset properties is via \code{\link{boxplotGenderGenotype}}, 
    \code{\link{boxplotGenderGenotypeBatch}} and \code{\link{scatterplotGenotypeWeight}} functions.
}
\section{Methods}{
    \describe{
        \item{testDataset}{\code{signature(object = "PhenList", equation = "character",depVariable = "character",
                pThreshold = "numeric",method = "character", keepList = "logical vector")}: Tests the fixed effects of the model}
        \item{scatterplotGenotypeWeight}{\code{signature(object = "PhenList", depVariable = "character")}: Creates a scatter plot body weight versus dependant variable}
        \item{boxplotGenderGenotypeBatch}{\code{signature(object = "PhenList", depVariable = "character")}: Creates a box plot split by gender, genotype and batch}
        \item{boxplotGenderGenotype}{\code{signature(object = "PhenList", depVariable = "character")}: Creates a box plot split by gender and genotype}              
    }
}
\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\link{PhenList}}}

\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"),
            testGenotype="Sparc/Sparc")
    class(test)
    test2 <- PhenList(dataset=read.csv("./PhenStat/extdata/MAAJ_CCdata.csv"),
            testGenotype="Arid4a/Arid4a",refGenotype="+/+",hemiGenotype="Arid4a/+",
            dataset.colname.weight="Weight.Value")
    test2$dataset.stat
}


\keyword{classes}