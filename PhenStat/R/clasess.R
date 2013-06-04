#    CLASSES.R defines classes that are used in PhenStat package: 
#    PhenList, PhenTestResult which extends LargeDataObject of limma package
require(methods)
setClass("LargeDataObject")
# LargeDataObject is the class from limma package developed by Gordon Smyth
# Here we are using LargeDataObject for the storage of phenotype data and for the storage of model fitting results
# and additional info 
printHead <- function(x)
#  Print leading 5 elements or rows of atomic object
#  Gordon Smyth
#  May 2003.  Last modified 14 April 2009.
{
    if(is.atomic(x)) {
        d <- dim(x)
        if(length(d)<2) which <- "OneD"
        if(length(d)==2) which <- "TwoD"
        if(length(d)>2) which <- "Array"
    } else {
        if(inherits(x,"data.frame")) {
            d <- dim(x)
            which <- "TwoD"
        } else {
            if(is.call(x))
            which <- "Call"
            else {
                if(is.recursive(x))
                which <- "Recursive"
                else
                which <- "Other"
            }
        }
    }
    switch(which,
            OneD={
                n <- length(x)
                if(n > 20) {
                    print(x[1:5])
                    cat(n-5,"more elements ...\n")
                } else
                print(x)
            },
            TwoD={
                n <- d[1]
                if(n > 10) {
                    print(x[1:5,])
                    cat(n-5,"more rows ...\n")
                } else
                print(x)
            },
            Array={
                n <- d[1]
                if(n > 10) {
                    dn <- dimnames(x)
                    dim(x) <- c(d[1],prod(d[-1]))
                    x <- x[1:5,]
                    dim(x) <- c(5,d[-1])
                    if(!is.null(dn[[1]])) dn[[1]] <- dn[[1]][1:5]
                    dimnames(x) <- dn
                    print(x)
                    cat(n-5,"more rows ...\n")
                } else
                print(x)
            },
            Recursive={
                n <- length(x)
                if(n) {
                    i <- names(x)
                    if(is.null(i)) i <- seq_len(n)
                    for (what in i) {
                        y <- x[[what]]
                        cat("$",what,"\n",sep="")
                        Recall(y)
                        cat("\n")
                    }
                }
            },
            Call=,Other=print(x)
            )
}

setMethod("show","LargeDataObject",
        #  Print and show method large data objects
        #  Gordon Smyth
        #  May 2003
        function(object) {
            cat("An object of class \"",class(object),"\"\n",sep="")
            for (what in names(object)) {
                x <- object[[what]]
                cat("$",what,"\n",sep="")
                printHead(x)
                cat("\n")
            }
            for (what in setdiff(slotNames(object),".Data")) {
                x <- slot(object,what)
                if(length(x) > 0) {
                    cat("@",what,"\n",sep="")
                    printHead(x)
                    cat("\n")
                }
            }
        })

setClass("PhenList",
        #  Linear model fit
        representation("list")
        )


setIs("PhenList","LargeDataObject")

dim.PhenList <- function(x) if(is.null(x$phendata)) c(0, 0) else dim(as.matrix(x$phendata))


length.PhenList <- function(x) prod(dim(x))

dimnames.PhenList <- function(x) dimnames(x$phendata)
assign("dimnames<-.PhenList",function(x,value)
        {
            dimnames(x$phendata) <- value
            x
        })

as.matrix.PhenList <- function(x,...) as.matrix(x$phendata)

setClass("PhenTestResult",
        #  model tests results
        representation("list")
        )

setIs("PhenTestResult","LargeDataObject")

dim.PhenTestResult <- function(x) if(is.null(x$modelOutput)) c(0, 0) else dim(as.matrix(x$modelOutput))


length.PhenTestResult <- function(x) prod(dim(x))

dimnames.PhenTestResult <- function(x) dimnames(x$modelOutput)
assign("dimnames<-.PhenTestResult",function(x,value)
        {
            dimnames(x$modelOutput) <- value
            x
        })

as.matrix.PhenTestResult <- function(x,...) as.matrix(x$modelOutput)

