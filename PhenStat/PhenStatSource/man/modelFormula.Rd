\name{modelFormula}
\alias{modelFormula}
\title{Method "modelFormula"}
\description{
    Creates formula for the start Mixed Model based on equation and number of sexes in the data.
}
\usage{    
    modelFormula(equation, numberofsexes, depVariable)
}
\arguments{
    \item{equation}{a character string defining the equation to use (withWeight or withoutWeight); mandatory argument}
    \item{numberofsexes}{a numerical value for levels of sex in the dataset (1 or 2); mandatory argument} 
    \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
}
\value{
    Returns formula to use.
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    formula <- modelFormula(equation="withWeight",
            numberofsexes=2, 
            depVariable="Lean.Mass")
}