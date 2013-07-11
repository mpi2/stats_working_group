\name{PhenList-class}
\docType{class}

\alias{PhenList-class}
\alias{testDataset,PhenList,character,character,numeric,character-method}
\alias{buildStartModel,PhenList,character,character,numerical,logical vector-method}
\alias{scatterplotGenotypeWeight,PhenList,character-method}
\alias{boxplotGenderGenotypeBatch,PhenList,character-method}
\alias{boxplotGenderGenotype,PhenList,character-method}  

\title{Class "PhenList"}

\description{A list-based S4 class for storing phenotypic data and support basic operations: quality check and statistics calculation}

\section{Create PhenList object}{
    The appropriate way to create PhenList object is via \code{\link{PhenList}} function.
}
\section{Dataset graphs}{
    The appropriate way to check the dataset properties is vie \code{\link{boxplotGenderGenotype}}, 
    \code{\link{boxplotGenderGenotypeBatch}} and \code{\link{scatterplotGenotypeWeight}} functions.
}
\section{Methods}{
    \describe{
        \item{testDataset}{\code{signature(object = "PhenList", equation = "character",depVariable = "character",
                pThreshold = "numeric",method = "character")}: Tests the fixed effects of the model}
        \item{buildStartModel}{\code{signature(object = "PhenList", equation = "character",depVariable = "character",
                        pThreshold = "numeric",keepList = "logical vector")}: Builds start model based in provided values 
            of significances for fixed effects}
        \item{scatterplotGenotypeWeight}{\code{signature(object = "PhenList", depVariable = "character")}: Creates a scatter plot body weight versus dependant variable}
        \item{boxplotGenderGenotypeBatch}{\code{signature(object = "PhenList", depVariable = "character")}: Creates a box plot split by gender, genotype and batch}
        \item{boxplotGenderGenotype}{\code{signature(object = "PhenList", depVariable = "character")}: Creates a box plot split by gender and genotype}
               
    }
}
\references{
   Placeholder for references
}
\author{Natalja Kurbatova}

\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    test2 <- PhenList(dataset=read.csv("./PhenStat/extdata/MAAJ_CCdata.csv"),dataset.colname.weight="Weight.Value",
            dataset.clean=TRUE,testGenotype="Arid4a/Arid4a",refGenotype="+/+",hemGenotype="Arid4a/+")
}


\keyword{classes}