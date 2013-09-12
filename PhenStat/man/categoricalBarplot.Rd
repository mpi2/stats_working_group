\name{categoricalBarplot}
\alias{categoricalBarplot}
\title{Method "categoricalBarplot"}
\description{
    Creates stacked bar plot(s) to compare proportions seen in a categorical variable between different genotypes. 
    Graphs are created for all data and also for each gender individually.
}
\usage{
    categoricalBarplot(phenList, phenTestResult, depVariable, graphingName=NULL)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; mandatory argument}
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class that comes from the function \code{testDataset}; mandatory argument}
    \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
    \item{graphingName}{a character string defining the name to be used on the graph for the dependent variable; defaults to NULL}
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenList}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/categorical_example1_aff3_Xray.csv"),
            testGenotype="Aff3/Aff3")
    result <- testDataset(test,depVariable="Thoracic.Processes")       
    categoricalBarplot(test, result, depVariable="Thoracic.Processes", graphingName="Thoracic Processes")
    
    }