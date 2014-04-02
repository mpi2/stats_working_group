\name{columnChecks}
\alias{columnChecks}
\title{Method "columnChecks"}
\description{
    Checks the particular column of the dataset for eligibility.
    Returns: presence of column, all data are numeric, number of levels that passed check 
    (number of data points for each genotype/sex combination is at least equals to threshold).
}
\usage{
    columnChecks(dataset, columnName, dataPointsThreshold=4)
}
\arguments{
    \item{dataset}{data frame created from file or from another source; mandatory argument}
    \item{columnName}{name of the column to check; mandatory argument}
    \item{dataPointsThreshold}{number of data points for each genotype/sex combination; default value is set to 4}
}

\value{
    Boolean vector of length 3: presence of column, all data are numeric, number of levels that passed check 
    (number of data points for each genotype/sex combination is at least equals to threshold). 
}
\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\keyword{internal}

