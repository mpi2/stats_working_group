\name{columnLevels}
\alias{columnLevels}
\title{Method "columnLevels"}
\description{
    Method removes null values from the column and calculates column basic stat.
    Returns: No of data points, No of levels, No of Genotype/Sex combinations, 
    No of data points for each one combination.
}
\usage{
    columnLevels(dataset, columnName)
}
\arguments{
    \item{dataset}{data frame created from file or from another source; 
    mandatory argument}
    \item{columnName}{name of the column to check; mandatory argument}
}

\value{Numeric vector of length 6 or 7: No of data points, No of levels, No of 
    Genotype/Sex combinations, 
    No of data points for each one combination.
}
\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and 
    Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} 
    \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide 
    using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\keyword{internal}

