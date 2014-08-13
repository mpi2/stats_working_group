\name{getName}
\alias{getName}
\title{Method "getName"}
\description{
    Returns IMPC object name by its ID. An internal function, if called directly please use with caution.
}
\usage{
    getName(fieldNameFrom,fieldNameTo,fieldValueFrom)
}
\arguments{
    \item{fieldNameFrom}{name of the IMPC field for the object ID; mandatory argument}
    \item{fieldNameTo}{name of the IMPC field for the object name; mandatory argument}
    \item{fieldValueFrom}{actual IMPC object's ID to get the name for; mandatory argument}
}
\value{
    Returns name of the object.
}

\author{Natalja Kurbatova, Jeremy Mason}
\examples{
   getName("pipeline_stable_id","pipeline_name","MGP_001")
}

