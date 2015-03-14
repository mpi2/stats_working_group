\name{determiningLambda}
\alias{determiningLambda}
\title{Method "determiningLambda"}
\description{
    Function role
    1.  Rescale data if needed (Box-Cox transformation can only be applied on data >0)
    2.  Calculate lambda value for a Box-Cox transformation 
    3.  Assess whether transformation is required (if confidence interval includes 1 then no transformation required)
    4.  Return output to allow application of the transformation if necessary using 
        the returned lambda value and any associated rescaling that was needed.
    5.  Requires a function that will be fitted determined from formulaAssessingBoxCox
}
\usage{
    determiningLambda(phenList, depVariable, equation="withWeight")
}
\arguments{
        \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; mandatory argument}       
        \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
        \item{equation}{a character string defining the equation to use. 
            Possible values "withWeight" (default),"withoutWeight"}
}
\value{
    Returns a vector with five elements:
        1. 95% confidence interval for lambda  
        2. Midepoint of the 95% confidence interval
        3. Calculated lambda value
        4. TransformationRequired: TRUE/FALSE
        5. Calculated scale shift value
}
\author{Natasha Karp}
\examples{
        file <- system.file("extdata", "test1.csv", package="PhenStat")
        test <- PhenList(dataset=read.csv(file),
                testGenotype="Sparc/Sparc")
        determiningLambda(test,"Lean.Mass",equation="withoutWeight")
}    
\keyword{internal}