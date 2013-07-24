\name{PhenList}
\alias{PhenList}
\title{Method "PhenList"}
\description{
    Creates PhenList object, checks the data integrity, rename colums when requested, 
    provide simple statistics ablut dataset. 
    Returns an instance of the PhenList object created from the data file
}
\usage{
    PhenList(dataset, outputMessages=TRUE, dataset.clean=TRUE, refGenotype='+/+', 
            testGenotype=NULL, hemiGenotype=NULL, dataset.colname.batch=NULL, dataset.colname.genotype=NULL, 
            dataset.colname.gender=NULL, dataset.colname.weight=NULL, dataset.values.missingValue=NULL, 
            dataset.values.male=NULL, dataset.values.female=NULL, dataset.stat=NULL)
    
}
\arguments{
    \item{dataset}{data frame created from file or from another source}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages}
    \item{dataset.clean}{flag: "FALSE" value for no modification of the dataset; "TRUE" value to clean dataset if needed}
    \item{testGenotype}{test genotype value}
    \item{hemiGenotype}{genotype value for heterozygos that will be changed to testGenotype value}
    \item{refGenotype}{reference genotype, assigned default value is "+/+"}    
    \item{dataset.colname.batch}{column name within "dataset" dataset for the batch effect}
    \item{dataset.colname.genotype}{column name within "dataset" dataset for the genotype info}
    \item{dataset.colname.gender}{column name within "dataset" dataset for the gender info}
    \item{dataset.colname.weight}{column name within "dataset" dataset for the weight info}
    \item{dataset.values.missingValue}{value used as missing value in the "dataset" dataset}
    \item{dataset.values.male}{value used to label "males" in the "dataset" dataset}
    \item{dataset.values.female}{value used to label "females" in the "dataset" dataset}
    \item{dataset.stat}{simple statistics obtained from "dataset" dataset}   
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
            dataset.clean=TRUE,testGenotype="Arid4a/Arid4a",refGenotype="+/+",hemiGenotype="Arid4a/+")
}