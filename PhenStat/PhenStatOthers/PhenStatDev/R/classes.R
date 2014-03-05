## Copyright © 2011-2013 EMBL - European Bioinformatics Institute
## 
## Licensed under the Apache License, Version 2.0 (the "License"); 
## you may not use this file except in compliance with the License.  
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##------------------------------------------------------------------------------
## CLASSES.R defines classes that are used in PhenStat package: 
## PhenList, PhenTestResult objects extends LargeDataObject of limma package
##------------------------------------------------------------------------------
setClass("PhenList",
        ##  Linear model fit
        representation("list")
        )

setIs("PhenList","LargeDataObject")

dim.PhenList <- function(x) 
if(is.null(x$dataset)) c(0, 0) else dim(as.matrix(x$dataset))

dimnames.PhenList <- function(x) dimnames(x$dataset)

assign("dimnames<-.PhenList",function(x,value)
        {
            dimnames(x$dataset) <- value
            x
        })

as.matrix.PhenList <- function(x,...) as.matrix(x$dataset)

#-----------------------------------------------------------------------------------
setClass("PhenTestResult",
        ##  model tests results
        representation("list")
        )

setIs("PhenTestResult","LargeDataObject")

method.PhenTestResult <- function(x) 
if(is.null(x$model.output)) NA else class(x$model.output)
##------------------------------------------------------------------------------