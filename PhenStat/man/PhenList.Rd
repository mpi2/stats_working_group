\name{PhenList}
\alias{PhenList}
\title{Method "PhenList"}
\description{
    Creates PhenList object, checks the data integrity, rename colums when requested, 
    provide simple statistics ablut dataset. 
    Returns an instance of the PhenList object created from the data file
}
\usage{
    PhenList(dataset=matrix(0,0,0), interactionMode=TRUE, refGenotype='+/+', dataset.stat=NULL, 
            dataset.colname.batch=NULL,dataset.colname.genotype=NULL, dataset.colname.gender=NULL, 
            dataset.colname.weight=NULL, dataset.values.missingValue=NULL,dataset.values.male=NULL, 
            dataset.values.female=NULL, dataset.clean=FALSE, testGenotype=NULL, hemGenotype=NULL) 
    
}
\arguments{
    \item{dataset}{data frame created from file or from another source}
    \item{interactionMode}{flag: "FALSE" value for the no interaction messages to the user; "TRUE" value to interact}
    \item{refGenotype}{reference genotype, assigned default value is "+/+"}
    \item{dataset.stat}{simple statistics obtained from "dataset" dataset}
    \item{dataset.colname.batch}{column name within "dataset" dataset for the batch effect}
    \item{dataset.colname.genotype}{column name within "dataset" dataset for the genotype info}
    \item{dataset.colname.gender}{column name within "dataset" dataset for the gender info}
    \item{dataset.colname.weight}{column name within "dataset" dataset for the weight info}
    \item{dataset.values.missingValue}{value used as missing value in the "dataset" dataset}
    \item{dataset.values.male}{value used to label "males" in the "dataset" dataset}
    \item{dataset.values.female}{value used to label "females" in the "dataset" dataset}
    \item{dataset.clean}{flag: "FALSE" value for no modification of the dataset; "TRUE" value to clean dataset if needed}
    \item{testGenotype}{test genotype}
    \item{hemGenotype}{genotype value for heterozygos that will be changed to testGenotype value}
}
\value{
    Returns an instance of the \code{\linkS4class{PhenList}} class.
}


\references{
    Publication placeholder
    
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}

\seealso{\code{\linkS4class{PhenTestResult}}}

\examples{
    test <- PhenList(dataset=read.csv("./PhenStat/extdata/MBAU_data.csv"))
    test2 <- PhenList(dataset=read.csv("./PhenStat/extdata/MAAJ_CCdata.csv"),dataset.colname.weight="Weight.Value",
            dataset.clean=TRUE,testGenotype="Arid4a/Arid4a",refGenotype="+/+",hemGenotype="Arid4a/+")
}