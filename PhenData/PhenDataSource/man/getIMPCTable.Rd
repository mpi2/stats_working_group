\name{getIMPCTable}
\alias{getIMPCTable}
\title{Method "getIMPCTable"}
\description{
    Returns the IMPC dataset. There are optional arguments to subset the IMPC dataset. 
    Be aware that without filters this file output will be very large and the query will be time intensive.
}
\usage{
    getIMPCTable(fileName="PhenData_IMPC",PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, ParameterID=NULL, 
            AlleleID=NULL, StrainID=NULL, multipleFiles=TRUE,recordsPerFile=10000)
}
\arguments{
    \item{fileName}{name of the file where to save resulting table with IMPC objects; mandatory argument}
    \item{PhenCenterName}{IMPC phenotyping center; mandatory argument}
    \item{PipelineID}{IMPC pipeline ID; mandatory argument}
    \item{ProcedureID}{IMPC procedure ID; mandatory argument}
    \item{ParameterID}{IMPC parameter ID; mandatory argument}
    \item{AlleleID}{IMPC allele ID; mandatory argument}
    \item{StrainID}{IMPC strain ID; optional argument}
    \item{multipleFiles}{flag: "FALSE" value to get all records into one specified file; 
        "TRUE" value (default) to split records across multiple files named starting with 'fileName' and 
        with digit at the end}
    \item{recordsPerFile}{number that specifies how many records to write into one file; default value is 10000}
}
\value{
    Returns the IMPC dataset.
}

\author{Natalja Kurbatova, Jeremy Mason}
\examples{
    getIMPCTable("./IMPCData_WTSI","WTSI","MGP_001","IMPC_CBC_001")
}
