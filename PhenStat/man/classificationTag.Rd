\name{classificationTag}
\alias{classificationTag}
\title{Method "classificationTag"}
\description{
    Returns classification tag for sexual dimorphism
}
\usage{
    classificationTag(result, interactionMode = TRUE, phenotypeThreshold = 0.01)
}
\arguments{
    \item{result}{instance of the \code{\linkS4class{PhenTestResult}} class}
    \item{interactionMode}{flag: "FALSE" value for the no interaction messages to the user; "TRUE" value to interact}
    \item{phenotypeThreshold}{a numerical value defining the threshold to use in classification tag}
}
\value{
    Returns model fitting results in a text format
}

\references{
    Placeholder for references
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    result <- testDataset(test,depVariable="Bone.Area")
    result <- buildFinalModel(test,result)
    classificationTag(result)
}