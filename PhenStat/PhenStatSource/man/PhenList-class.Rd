\name{PhenList-class}
\docType{class}
\alias{PhenList-class}
\alias{dataset,PhenList-method}
\alias{refGenotype,PhenList-method}
\alias{testGenotype,PhenList-method}
\alias{hemiGenotype,PhenList-method}
\alias{getColumn,PhenList,character-method}
\alias{getColumnBatchAdjusted,PhenList,character-method}
\alias{getStat,PhenList-method}
\alias{getVariables,PhenList-method}
\alias{multipleBatches,PhenList-method}
\alias{setBatch,PhenList,character-method}
\alias{setGenotype,PhenList,character-method}
\alias{setMissingValue,PhenList,character-method}
\alias{setSex,PhenList,character-method}
\alias{setWeight,PhenList,character-method}
\title{Class "PhenList"}
\description{A list-based S4 class for storing phenotypic data. Helps to support data integrity checks and statistics calculation.
The PhenList object can be created by using function \code{\link{PhenList}}.}

\section{Explore PhenList object}{
    PhenList object instance contains the following slots:
    
    1. "datasetPL" where cleaned and checked dataset is stored: \code{getDataset(phenList)}
    
    2. "testGenotype" where the genotype level to test is stored:  \code{testGenotype(phenList)}
    
    3. "refGenotype" where reference genotype  value is stored with default value set to "+/+":  
    
    \code{refGenotype(phenList)}
    
    4. "hemiGenotype" if defined contains the genotype value for hemizygous: \code{hemiGenotype(phenList)}
    
    6. "dataset.clean" flag value is stored,  see \code{\link{PhenList}} for more details:
    
    \code{phenList@dataset.clean}
    
    7. "dataset.colname" if defined contains dataset column names that have been renamed:
    
        - \code{phenList@dataset.colname.batch} column name for batch values
       
        - \code{phenList@dataset.colname.genotype} column name for genotype values
    
        - \code{phenList@dataset.colname.sex} column name for sex values
    
        - \code{phenList@dataset.colname.weight} column name for weight values
        
    8. "dataset.values" if defined contains dataset values that have been changed during dataset cleaning process:
    
        - \code{phenList@dataset.values.missingValue} value used as missing value in the orginal dataset
        
        - \code{phenList@dataset.values.male} value used to label "males" in the original dataset
        
        - \code{phenList@dataset.values.female} value used to label "females" in the original dataset
}
\section{Slots}{
    \describe{
        \item{\code{datasetPL}:}{Object of class \code{"data.frame"} ~~ dataset to work with}
        \item{\code{refGenotype}:}{Object of class \code{"character"} ~~ reference genotype}
        \item{\code{testGenotype}:}{Object of class \code{"character"} ~~ test genotype}
        \item{\code{hemiGenotype}:}{Object of class \code{"character"} ~~ hemi genotype}
        \item{\code{dataset.colname.batch}:}{Object of class \code{"character"} ~~ column name for batch values}
        \item{\code{dataset.colname.genotype}:}{Object of class \code{"character"} ~~ column name for genotype values}
        \item{\code{dataset.colname.sex}:}{Object of class \code{"character"} ~~ column name for sex values}
        \item{\code{dataset.colname.weight}:}{Object of class \code{"character"} ~~ column name for weigth values}
        \item{\code{dataset.values.missingValue}:}{Object of class \code{"character"} ~~ value used as missing value 
            in the orginal dataset}
        \item{\code{dataset.values.male}:}{Object of class \code{"character"} ~~ value used to label "males" 
            in the original dataset}
        \item{\code{dataset.values.female}:}{Object of class \code{"character"} ~~ value used to label "females" 
            in the original dataset}
        \item{\code{dataset.clean}:}{Object of class \code{"logical"} ~~ flag value is stored}
    }
}
\section{Methods}{
    \describe{
        \item{getDataset}{(accessor): Returns dataset}
        \item{refGenotype}{(accessor): Returns reference genotype }
        \item{testGenotype}{(accessor): Returns test genotype }
        \item{hemiGenotype}{(accessor): Returns hemi genotype if specified}
        \item{getColumn}{Returns specified column if exists}
        \item{getColumnBatchAdjusted}{Returns specified column adjusted for batch effect}
        \item{getStat}{Returns simple dataset statistics}
        \item{getVariables}{Returns dataset column names}
        \item{batchIn}{Returns TRUE if the batch is in the dataset, FALSE otherwise}
        \item{weightIn}{Returns TRUE if the weight is in the dataset, FALSE otherwise}
        \item{multipleBatches}{Returns TRUE if the batches are variable in the dataset, FALSE otherwise}
        \item{noSexes}{Returns number of sexes in the dataset}
        \item{setBatch}{Sets the batch column - change the column names}
        \item{setGenotype}{Sets the genotype column - change the column names}
        \item{setMissingValue}{Sets the missing value}
        \item{setSex}{Sets the sex column - change the column names}
        \item{setWeight}{Sets the weight column - change the column names}
    }
}

\author{Natalja Kurbatova}

\seealso{\code{\link{PhenList}}}

\examples{
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenList(dataset=read.csv(file),
            testGenotype="Sparc/Sparc")
    class(test)
    
    file <- system.file("extdata", "test2.csv", package="PhenStat")
    test2 <- PhenList(dataset=read.csv(file),
            testGenotype="Arid4a/Arid4a",
            refGenotype="+/+",
            hemiGenotype="Arid4a/+",
            dataset.colname.weight="Weight.Value")
    getStat(test2)
    testGenotype(test2)
    refGenotype(test2)
    hemiGenotype(test2)
}

\keyword{methods}
\keyword{classes}