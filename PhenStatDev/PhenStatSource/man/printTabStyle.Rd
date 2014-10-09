\name{printTabStyle}
\alias{printTabStyle}
\title{Method "printTabStyle"}
\description{
    Makes ready list of strings to be printed in tabular form, creates a 
    table row from it. 
    Adds empty positions if needed in order to make all strings the same by 
    length and adds
    "tabSep" character between strings.
    Returns text that is ready to be printed (one row of a table). 
}
\usage{
    printTabStyle(textList,positions,tabSep="|")
}
\arguments{
    \item{textList}{list of strings to prepare}
    \item{positions}{desired length of one cell in the table}
    \item{tabSep}{character used to separate columns in the table; 
    default value is "|"}
}
\value{Returns text that is ready to be printed (one row of a table). 
}
\author{Natalja Kurbatova}
\examples{
    row_sep <- rep("----------",3)
    print(printTabStyle(row_sep,11))
    list1 <- c("Owner","Model","Cc")
    print(printTabStyle(list1,11))
    print(printTabStyle(row_sep,11))
    list2 <- c("Natalja","Honda","125")
    print(printTabStyle(list2,11))
    print(printTabStyle(row_sep,11))
}    
\keyword{internal}