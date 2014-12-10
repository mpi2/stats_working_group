\name{recommendMethod}
\alias{recommendMethod}
\title{Method "recommendMethod"}
\description{     
   Checks the dependent variable and dataset to make choose the appropriate statistical method. 
   Returns the method or the list of methods that are appropriate for statistical analysis. 
}    
\usage{
    recommendMethod(phenList=NULL, depVariable=NULL, outputMessages=TRUE)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; mandatory argument}       
    \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages; 
        default value TRUE}        
}
\value{
    Returns the method or the list of methods that are appropriate for statistical analysis.
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\seealso{\code{\linkS4class{PhenList}}}

\examples{
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    recommendMethod(test,"Lean.Mass")
}