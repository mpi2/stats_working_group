\name{PhenList}
\alias{PhenList}
\title{Returns an instance of the PhenList object created from the data file.}
\description{
    Creates PhenList object, checks the data integrity, rename colums when requested, provide simple statistics ablut dataset. 
}
\usage{
    PhenList(dataFrame)
}
\arguments{
    \item{dataFrame}{data frame created from file or from another source}
}
\value{
    Returns an instance of the \code{\linkS4class{PhenList}} class.
}
\references{
    Publication placeholder

}
\author{Natalja Kurbatova}

\seealso{\code{\linkS4class{PhenList}}}

\examples{
    test <- PhenList(phendata=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    test2 <- PhenList(phendata=read.csv("./PhenStat/extdata/MAAJ_CCdata.csv"),weight.colname="Weight.Value")
}