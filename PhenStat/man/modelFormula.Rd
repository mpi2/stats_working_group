\name{modelFormula}
\alias{modelFormula}
\title{Method "modelFormula"}
\description{
    Create formula for the start Mixed Model based on equation and number of genders in the data
}
\usage{    
    modelFormula(equation, numberofgenders, depVariable)
}
\arguments{
    \item{equation}{a character string defining the equation to use}
    \item{numberofgenders}{a numerical value for levels of gender in the dataset} 
    \item{depVariable}{a character string defining the dependent variable of interest}
}
\value{
    Return start model
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    formula <- modelFormula(equation="withWeight",numberofgenders=2, 
            depVariable="Bone.Area")
}