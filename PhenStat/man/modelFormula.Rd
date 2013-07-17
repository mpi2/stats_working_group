\name{modelFormula}
\alias{modelFormula}
\title{Method "modelFormula"}
\description{
    Creates formula for the start model based on equation and number of genders in the data
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
    Returns start model
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    formula <- modelFormula(equation="withWeight",numberofgenders=2, depVariable="Bone.Area")
}