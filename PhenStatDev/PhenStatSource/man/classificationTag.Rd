\name{classificationTag}
\alias{classificationTag}
\title{Method "classificationTag"}
\description{
    Returns a classification tag to assign a sexual dimorphism assessment of 
    the phenotypic change.
}
\usage{
    classificationTag(phenTestResult, userMode = "summaryOutput", 
            phenotypeThreshold = 0.01,outputMessages=TRUE)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} 
    class that comes from the function 
        \code{\link{testDataset}}; mandatory argument}
    \item{userMode}{flag: "vectorOuput" a sexual dimorphic classification is 
    assigned with a proviso of later assessing 
        the overall statistical significance; "summaryOutput" the 
        phenotypeThreshold is used to assess the overall statistical 
        significance and then if significant the sexual dimorphic 
        classification determined; defaults to summaryOutput}
    \item{phenotypeThreshold}{a numerical value defining the threshold to 
    use in classificationTag in determining 
        whether the genotype effect is classed as significant or not; 
        default value 0.01}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; 
    "TRUE" value to show output messages
        ; default value TRUE}        
}
\value{Returns a classification tag to assign a sexual dimorphism assessment 
of the phenotypic change. In the case of Fisher Exact Test it is not sexual 
dimorphism classification, but rather the overall estimation of the signals 
significance across the different datasets.   
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and 
    Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} 
    \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410
    
    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide 
    using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    result <- testDataset(test,
            depVariable="Lean.Mass")
    classificationTag(result, 
            userMode="summaryOutput", 
            phenotypeThreshold=0.001)
    classificationTag(result, 
            userMode="vectorOutput")
}