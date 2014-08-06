\name{getName}
\alias{getName}
\title{Method "getName"}
\description{
    Returns IMPC object name by its ID.
}
\usage{
    getName(fieldNameFrom,fieldNameTo,fieldValueFrom)
}
\arguments{
    \item{fieldNameFrom}{name of the IMPC field for the object ID, for instance "pipeline_stable_id"; mandatory argument}
    \item{fieldNameTo}{name of the IMPC field for the object name, for instance "pipeline_name"; mandatory argument}
    \item{fieldValueFrom}{actual IMPC object''s ID to get the name for, for instance "MGP_001"; mandatory argument}
}
\value{
    Returns name of the object.
}

\author{Natalja Kurbatova, Jeremy Mason}
\example{
   getName("pipeline_stable_id","pipeline_name","MGP_001")
}