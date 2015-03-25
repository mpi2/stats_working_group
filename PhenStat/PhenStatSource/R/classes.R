## Copyright © 2012-2014 EMBL - European Bioinformatics Institute
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
## PhenList, PhenTestResult objects
##------------------------------------------------------------------------------
setClass("PhenList",
        ##  Linear model fit
        representation(datasetPL = "data.frame",
                refGenotype = "character",
                testGenotype = "character",
                hemiGenotype = "character",
                dataset.colname.batch = "character",
                dataset.colname.genotype = "character",
                dataset.colname.sex = "character",
                dataset.colname.weight = "character",
                dataset.values.missingValue = "character",
                dataset.values.male = "character",
                dataset.values.female = "character",
                dataset.clean = "logical")
)
##------------------------------------------------
# Dimension functions
dim.PhenList <- function(x){ 
    if(is.null(x@datasetPL)) c(0, 0) 
    else 
    c(nrow(x@datasetPL),ncol(x@datasetPL))
}
#dimnames.PhenList <- function(x) dimnames(x@dataset)
##------------------------------------------------
# Accessor functions
getDataset = function(obj) obj@datasetPL
refGenotype = function(obj) obj@refGenotype
testGenotype = function(obj) obj@testGenotype
hemiGenotype = function(obj) obj@hemiGenotype
##------------------------------------------------
# Show function
setMethod(show, signature(object = "PhenList"),
        function(object){
            print(paste("Dataset dimension - ",dim(object)[1]," x ",dim(object)[2],sep=""))
            print(paste("Reference Genotype - ",refGenotype(object),sep=""))
            print(paste("Test Genotype - ",testGenotype(object),sep=""))
            #print(paste("Original Male/Female values - ",
            #ifelse(is.na(obj@dataset.values.male),"Male",obj@dataset.values.male),
            #"/",
            #ifelse(is.na(obj@dataset.values.female),"Female",obj@dataset.values.female) 
            #,sep=""))
            getStat(object)
        }
)
##------------------------------------------------
# Number of sexes
setGeneric("noSexes",
        function(obj)
        standardGeneric("noSexes"))
setMethod("noSexes", signature(obj = "PhenList"),
        function(obj){
            length(levels(getDataset(obj)$Sex))
        }
)
##------------------------------------------------
# Multiple Batches - TRUE or FALSE
setGeneric("multipleBatches",
        function(obj)
        standardGeneric("multipleBatches"))
setMethod("multipleBatches", signature(obj = "PhenList"),
        function(obj){
            if ('Batch' %in% colnames(getDataset(obj))){
                batchColumn <- na.omit(getColumn(obj,"Batch"))
                if (length(levels(batchColumn)) > 1)
                    TRUE
                else
                    FALSE
            }
            else 
                FALSE
        }
)
##------------------------------------------------
# Weight is in the dataset - TRUE or FALSE
setGeneric("weightIn",
        function(obj)
        standardGeneric("weightIn"))
setMethod("weightIn", signature(obj = "PhenList"),
        function(obj){
            if ('Weight' %in% colnames(getDataset(obj)))
                TRUE
            else 
                FALSE
        }
)
##------------------------------------------------
# Batch is in the dataset - TRUE or FALSE
setGeneric("batchIn",
        function(obj)
        standardGeneric("batchIn"))
setMethod("batchIn", signature(obj = "PhenList"),
        function(obj){
            if ('Batch' %in% colnames(getDataset(obj)))
                TRUE
            else 
                FALSE
        }
)
##------------------------------------------------
# Stat - new generic function
setGeneric("getStat",
        function(obj)
        standardGeneric("getStat"))

setMethod("getStat", signature(obj = "PhenList"),
        function(obj){
            ## Calculate statistics  
            dataset.stat <- data.frame(
                    Variables = colnames(getDataset(obj)),
                    Numeric = sapply(getDataset(obj), is.numeric),    
                    Continuous = sapply(getDataset(obj), function(x) if(is.numeric(x)) {
                                if (length(unique(x))/length(x)>0.05) TRUE else FALSE} else FALSE),
                    Levels = sapply(getDataset(obj), function(x) length(unique(x)) ),
                    NObs = sapply(getDataset(obj), function(x) length(na.omit(x))),
                    Mean = sapply(getDataset(obj), function(x) 
                            if(is.numeric(x) && (length(unique(x))/length(x)>0.05)) round(mean(na.omit(x)),digits=2) else NA),
                    StdDev = sapply(getDataset(obj), function(x) 
                            if(is.numeric(x) && (length(unique(x))/length(x)>0.05)) round(sd(na.omit(x)),digits=2) else NA),
                    Minimum = sapply(getDataset(obj), function(x) 
                            if(is.numeric(x) && (length(unique(x))/length(x)>0.05)) round(min(na.omit(x)),digits=2) else NA),
                    Maximum = sapply(getDataset(obj), function(x) 
                            if(is.numeric(x) && (length(unique(x))/length(x)>0.05)) round(max(na.omit(x)),digits=2) else NA)
                    )
            rownames(dataset.stat) <- NULL
            print(dataset.stat,quote = FALSE)
        }
)
##------------------------------------------------
# Column of interest
setGeneric("getColumn",
        function(obj,columnName)
        standardGeneric("getColumn"))
setMethod("getColumn", signature(obj = "PhenList",columnName="character"),
        function(obj, columnName){
            if (c(columnName) %in% colnames(getDataset(obj))){
                columnOfInterest <- getDataset(obj)[,c(columnName)]               
            }
            else {
                columnOfInterest <- NULL
            }
            columnOfInterest
        }
)
##------------------------------------------------
# Column of interest adjusted for batch
setGeneric("getColumnBatchAdjusted",
        function(obj,columnName)
        standardGeneric("getColumnBatchAdjusted"))
setMethod("getColumnBatchAdjusted", signature(obj = "PhenList",columnName="character"),
        function(obj, columnName){
            if (c(columnName) %in% colnames(getDataset(obj))){
                columnOfInterest <- getDataset(obj)[,c(columnName)] 
                if (multipleBatches(obj)){
                    # Adjusted for batch depVariable values WITHOUT transformation!
                    nullFormula=as.formula(paste(columnName, "~", paste("1", sep="+")))        
                    model.null=do.call("lme", args = list(nullFormula, random=~1|Batch, getDataset(obj), 
                                    na.action="na.exclude", method="ML"))
                    # Find the depVariable column adjusted for the batch effect
                    columnOfInterest=resid(model.null)     
                }         
            }
            else {
                columnOfInterest <- NULL
            }
            columnOfInterest
        }
        )
##------------------------------------------------
# Variables
setGeneric("getVariables",
        function(obj)
        standardGeneric("getVariables"))
setMethod("getVariables", signature(obj = "PhenList"),
        function(obj){
            colnames(getDataset(obj))
        }
        )
##------------------------------------------------
# set functions
setGeneric("setBatch",
        function(obj,dataset.colname.batch)
        standardGeneric("setBatch"))

setMethod("setBatch", signature(obj = "PhenList",dataset.colname.batch="character"),
        function(obj, dataset.colname.batch){
            datasetCopy <- getDataset(obj)
            if(dataset.colname.batch!='Batch'){
                # check here for the existing column named 'Batch'
                if ('Batch' %in% colnames(datasetCopy)){
                    colnames(datasetCopy)[colnames(datasetCopy) == 'Batch'] <-'Original.Batch'   
                }
                colnames(datasetCopy)[colnames(datasetCopy) == dataset.colname.batch] <-'Batch'
            }
            obj@datasetPL <- datasetCopy
            obj@dataset.colname.batch <- dataset.colname.batch
            obj
        }
)
##------------------------------------------------
setGeneric("setGenotype",
        function(obj,dataset.colname.genotype)
        standardGeneric("setGenotype"))

setMethod("setGenotype", signature(obj = "PhenList",dataset.colname.genotype="character"),
        function(obj, dataset.colname.genotype){
            datasetCopy <- getDataset(obj)
            if(dataset.colname.genotype!='Genotype'){
                # check here for the existing column named 'Genotype'
                if ('Genotype' %in% colnames(datasetCopy)){
                    colnames(datasetCopy)[colnames(datasetCopy) == 'Genotype'] <-'Original.Genotype'   
                }
                colnames(datasetCopy)[colnames(datasetCopy) == dataset.colname.genotype] <-'Genotype'
            }
            obj@datasetPL <- datasetCopy
            obj@dataset.colname.genotype <- dataset.colname.genotype
            obj
        }
)
##------------------------------------------------
setGeneric("setSex",
        function(obj,dataset.colname.sex)
        standardGeneric("setSex"))

setMethod("setSex", signature(obj = "PhenList",dataset.colname.sex="character"),
        function(obj, dataset.colname.sex){
            datasetCopy <- getDataset(obj)
            if(dataset.colname.sex!='Sex'){
                # check here for the existing column named 'Sex'
                if ('Sex' %in% colnames(datasetCopy)){
                    colnames(datasetCopy)[colnames(datasetCopy) == 'Sex'] <-'Original.Sex'   
                }
                colnames(datasetCopy)[colnames(datasetCopy) == dataset.colname.sex] <-'Sex'
            }
            obj@datasetPL <- datasetCopy
            obj@dataset.colname.sex <- dataset.colname.sex
            obj
        }
)
##------------------------------------------------
setGeneric("setWeight",
        function(obj,dataset.colname.weight)
        standardGeneric("setWeight"))

setMethod("setWeight", signature(obj = "PhenList",dataset.colname.weight="character"),
        function(obj, dataset.colname.weight){
            datasetCopy <- getDataset(obj)
            if(dataset.colname.weight!='Weight'){
                # check here for the existing column named 'Weight'
                if ('Sex' %in% colnames(datasetCopy)){
                    colnames(datasetCopy)[colnames(datasetCopy) == 'Weight'] <-'Original.Weight'   
                }
                colnames(datasetCopy)[colnames(datasetCopy) == dataset.colname.weight] <-'Weight'
            }
            obj@datasetPL <- datasetCopy
            obj@dataset.colname.weight <- dataset.colname.weight
            obj
        }
)
##------------------------------------------------
setGeneric("setMissingValue",
        function(obj,dataset.values.missingValue)
        standardGeneric("setMissingValue"))

setMethod("setMissingValue", signature(obj = "PhenList",dataset.values.missingValue="character"),
        function(obj, dataset.values.missingValue){
            datasetCopy <- getDataset(obj)
            ## Replace missing values specified in the user format with NA 
            datasetCopy[datasetCopy == dataset.values.missingValue] <- NA        
            obj@dataset <- datasetCopy
            if (is.na(obj@dataset.values.missingValue)){
                obj@dataset.values.missingValue <- dataset.values.missingValue
            }
            else{
                if (obj@dataset.values.missingValue!=dataset.values.missingValue)
                obj@dataset.values.missingValue <- c(obj@dataset.values.missingValue,dataset.values.missingValue)
            }
            obj
        }
)
#-----------------------------------------------------------------------------------
setClass("PhenTestResult",
        representation(
                analysedDataset = "data.frame",
                transformationRequired = "logical",
                lambdaValue = "numeric",
                scaleShift = "numeric",
                depVariable = "character",
                refGenotype = "character",
                testGenotype = "character",
                method = "character",
                parameters = "matrix",
                analysisResults="list")
        
)
# Stores statistical analysis results:  depending on method they can be gls, lme (MM and TF), loigstf (LR), 
# list of htestPhenStat objects (FE and RR) objects, plus common meta data. 
# 'analysedDataset' slot contains only columns of original dataset that have been used in the analysis
##------------------------------------------------
# Accessor functions
getVariable = function(obj) obj@depVariable
refGenotype = function(obj) obj@refGenotype
testGenotype = function(obj) obj@testGenotype
parameters = function(obj) obj@parameters
method = function(obj) obj@method
methodText = function(obj) 
        ifelse(method(obj) == "RR", "Reference Range Plus", 
        ifelse(method(obj) == "FE", "Fisher Exact Test", 
        ifelse(method(obj) == "MM", "Linear Mixed Model",
        ifelse(method(obj) == "LR", "Logistic Regression",
                        "Time as Fixed Effect"))))
analysisResults = function(obj) obj@analysisResults
analysedDataset = function(obj) obj@analysedDataset
transformationText = function(obj) 
        ifelse(obj@transformationRequired,
            paste(ifelse((obj@lambdaValue!=0),
                        paste(", power transformed with lambda value =",obj@lambdaValue)
                        ,", log transformed"),
                ifelse(!(obj@scaleShift==0),
                        paste(" and scale shift =",obj@scaleShift)
                        ,""))
        ,"")
transformation = function(obj) 
    ifelse(obj@transformationRequired,
        paste("lambda=",obj@lambdaValue,", scaleShift=",obj@scaleShift,sep="")
        ,"lambda=NA, scaleShift=NA")
##------------------------------------------------
# Number of sexes
setMethod("noSexes", signature(obj = "PhenTestResult"),
        function(obj){
            x <- analysedDataset(obj)
            length(levels(x$Sex))
        }
)
##------------------------------------------------
# Batch is in the dataset - TRUE or FALSE
setMethod("batchIn", signature(obj = "PhenTestResult"),
        function(obj){
            x <- analysedDataset(obj)
            if ('Batch' %in% colnames(x))
                TRUE
            else 
                FALSE
        }
)
##------------------------------------------------
# Weight is in the dataset - TRUE or FALSE
setMethod("weightIn", signature(obj = "PhenTestResult"),
        function(obj){
            x <- analysedDataset(obj)
            if ('Weight' %in% colnames(x))
                TRUE
            else 
                FALSE
        }
)
##------------------------------------------------
setGeneric("getCountMatrices",
        function(obj)
        standardGeneric("getCountMatrices"))

setMethod("getCountMatrices", signature(obj = "PhenTestResult"),
        function(obj){
            x <- analysedDataset(obj)
            noSexes <- length(levels(x$Sex))
            modeloutput <- analysisResults(obj)
            if (method(obj)=="RR"){
                for (i in seq_along(modeloutput)) {
                    val <- modeloutput[[i]]
                    if (comparison(val)=="High vs Normal/Low") {
                        if (analysedSubset(val)=="all"){
                            high_all_val <- matrixCount(val)["High",]
                            normallow_all_val <- matrixCount(val)["Normal/Low",]
                        }
                        if (analysedSubset(val)=="males"){
                            high_male_val <- matrixCount(val)["High",]    
                            normallow_male_val <- matrixCount(val)["Normal/Low",]                     
                        }
                        if (analysedSubset(val)=="females"){
                            high_female_val <- matrixCount(val)["High",] 
                            normallow_female_val <- matrixCount(val)["Normal/Low",]                          
                        }
                    }
                    else {
                        if (analysedSubset(val)=="all"){
                            low_all_val <- matrixCount(val)["Low",]     
                        }
                        if (analysedSubset(val)=="males"){
                            low_male_val <- matrixCount(val)["Low",]                         
                        }
                        if (analysedSubset(val)=="females"){
                            low_female_val <- matrixCount(val)["Low",]                           
                        }
                    }
                }
                all_matrix <- rbind(low_all_val,normallow_all_val-low_all_val,high_all_val)
                rownames(all_matrix) <- c("Low","Normal","High")
                if(noSexes==2){
                    male_matrix <- rbind(low_male_val,normallow_male_val-low_male_val,high_male_val)
                    female_matrix <- rbind(low_female_val,normallow_female_val-low_female_val,high_female_val) 
                    rownames(male_matrix) <- c("Low","Normal","High")
                    rownames(female_matrix) <- c("Low","Normal","High")
                    list(all=all_matrix,male=male_matrix,female=female_matrix)                   
                }
                else {
                    list(all=all_matrix)
                }
            }
            else if (method(obj)=="FE"){
                for (i in seq_along(modeloutput)) {
                    val <- modeloutput[[i]]
                    if (analysedSubset(val)=="all"){
                        all_matrix <- val@matrixCount  
                        rownames(all_matrix) <- rownames(val@matrixCount) 
                    }
                    if (analysedSubset(val)=="males"){
                        male_matrix <- val@matrixCount  
                        rownames(male_matrix) <- rownames(val@matrixCount)                    
                    }
                    if (analysedSubset(val)=="females"){
                        female_matrix <- val@matrixCount    
                        rownames(female_matrix) <- rownames(val@matrixCount)                        
                    }
                }
                if(noSexes==2){
                    list(all=all_matrix,male=male_matrix,female=female_matrix)                   
                }
                else {
                    list(all=all_matrix)
                }
                
            }
            else {
                NULL
            }        
        }
)
##------------------------------------------------
setGeneric("getGenotypeEffect",
        function(obj)
        standardGeneric("getGenotypeEffect"))

setMethod("getGenotypeEffect", signature(obj = "PhenTestResult"),
        function(obj){
            if (method(obj) %in% c("MM","TF","LR")){
                effect_values <- c(obj@analysisResults$model.output.summary["genotype_estimate"],
                        obj@analysisResults$model.output.summary["genotype_estimate_SE"])
                if (obj@transformationRequired)
                    effect_values <- reverseTransformValues(effect_values,obj@lambdaValue,obj@scaleShift)
                
                as.numeric(effect_values)
            }
            else {
                NULL
            }
        }
)
##------------------------------------------------
# Show function
setMethod(show, signature(object = "PhenTestResult"),
        function(object){
            x <- analysedDataset(object)
            noSexes <- length(levels(x$Sex))
            
            cat("****Information****\n")
            cat(paste("Reference Genotype: ", refGenotype(object),
                            "; Test Genotype: ", testGenotype(object),"\n",sep=""))
            cat(paste("Variable: ", getVariable(object),
                            transformationText(object),
                            "\n",sep=""))
           
            cat("\n******Method******\n")
            cat(paste(methodText(object),"\n",sep=""))
            
            if (length(parameters(object))>0){
                cat("\n****Parameters****\n")
                print(parameters(object),quote=FALSE)
            }
            if (method(object) %in% c("RR")){
                cat("\n*****Result*****\n")
                #Prepare data               
                nl <- data.frame(nr=c(1,2,3))
                for (i in seq_along(analysisResults(object))) {
                    val <- analysisResults(object)[[i]]
                    if (comparison(val)=="High vs Normal/Low"){
                        nl[paste("new",i)]<-as.vector(getColumnView(val))                  
                    }
                }
                nl<-nl[ , -which(names(nl) %in% c("nr"))]
                nh <- data.frame(nr=c(1,2,3))
                for (i in seq_along(analysisResults(object))) {
                    val <- analysisResults(object)[[i]]
                    if (comparison(val)=="Low vs Normal/High"){
                        nh[paste("new",i)]<-as.vector(getColumnView(val))                  
                    }
                }
                nh<-nh[ , -which(names(nh) %in% c("nr"))]
                if (noSexes==2){
                    colnames(nl) <- head(nl,1)
                    nl <- nl[-1,]
                    rownames(nl) <- c("p-value","ES")
                    cat("\n1) High vs Normal/Low\n")
                    print(nl,quote=FALSE)
                    colnames(nh) <- head(nh,1)
                    nh <- nh[-1,]
                    rownames(nh) <- c("p-value","ES")
                    cat("\n2) Low vs Normal/High\n")
                    print(nh,quote=FALSE)
                }
                else {
                    all <- cbind(nl,nh)
                    all <- all[-1,]
                    rownames(all) <- c("p-value","ES")
                    colnames(all) <- c("High vs Normal/Low","Low vs Normal/High")
                    print(all,quote=FALSE)
                }
                
            }
            if (method(object) %in% c("FE")){
                cat("\n*****Result*****\n")
                nl <- data.frame(nr=c(1,2,3))
                for (i in seq_along(analysisResults(object))) {
                    val <- analysisResults(object)[[i]]
                    nl[paste("new",i)]<-as.vector(getColumnView(val)) 
                }
                nl<-nl[ , -which(names(nl) %in% c("nr"))]
                if (noSexes==2){
                    colnames(nl) <- head(nl,1)
                    nl <- nl[-1,]
                    rownames(nl) <- c("p-value","Effect Size")
                    print(nl,quote=FALSE)
                }
                else {
                    cat(paste("p-value: ",nl[2],"\n",sep=""))
                    cat(paste("Effect Size: ",nl[3],"\n",sep=""))
                }
            }  
            if (method(object) %in% c("MM","TF","LR")) {       
                cat(paste("Equation: ",object@analysisResults$equation ,"\n",sep=""))
                        
                cat("\n*****Result*****\n")
                printLROutput(object)
            }
            #else 
            #  cat(paste("Method: ", method(object),
            #            "; Dataset: ", subsetText(object),sep=""))
            #cat("\n\nModel Output:\n")
            #write.table(result,quote=FALSE,col.names = FALSE)
            #cat("\nParameters:\n")
            #d<-lapply(test@parameters, function(x){cat(x); cat("\n")})
        }
)
##------------------------------------------------
# htest objects that come from FE and RR methods
setOldClass("htest")
setClass("htestPhenStat",
        representation(
                modelOutput = "htest",
                analysedSubset = "character",
                comparison = "character",
                ES = "numeric",
                matrixCount = "matrix")
)
##------------------------------------------------
# Accessor functions
pvalue = function(obj) obj@modelOutput$p.value
analysedSubset = function(obj) obj@analysedSubset
subsetText = function(obj) ifelse(analysedSubset(obj) == "all", "All", 
        ifelse(analysedSubset(obj) == "males", "Males only", "Females only"))
comparison = function(obj) obj@comparison
matrixCount = function(obj) obj@matrixCount
##------------------------------------------------
# Show method for htest 
setMethod(show, signature(object = "htestPhenStat"),
        function(object){
            if (length(comparison(object))>0){
                result <- matrix(0,4,1)
                result[1] <- subsetText(object)
                result[2] <- comparison(object)
                result[3] <- sprintf("%.4f",round(pvalue(object), digits=4))
                result[4] <-  paste(format(object@ES, nsmall = 0),"%",sep="")
                rownames(result) <- c("Analized subset:","Comparison:",
                        "p-value:","Effect size:")
            }
            else {
                result <- matrix(0,3,1)
                result[1] <- subsetText(object)
                result[2] <- sprintf("%.4f",round(pvalue(object), digits=4))
                result[3] <-  paste(format(object@ES, nsmall = 0),"%",sep="")
                rownames(result) <- c("Analised subset:","p-value:","Effect size:")
            }
            colnames(result) <- c("")
            
            print(result,quote=FALSE)
            #cat("Information:\n")
            #cat(paste("Variable: ", Variable(object),
            #          "; Reference Genotype: ", refGenotype(object),
            #          "; Test Genotype: ", testGenotype(object),"\n",sep=""))
            #if (method(object)=="RR")
            
            #else 
            #  cat(paste("Method: ", method(object),
            #            "; Dataset: ", subsetText(object),sep=""))
            #cat("\n\nModel Output:\n")
            #write.table(result,quote=FALSE,col.names = FALSE)
            #cat("\nParameters:\n")
            #d<-lapply(test@parameters, function(x){cat(x); cat("\n")})
        }
)
##------------------------------------------------
setGeneric("getMatrix",
        function(obj,phenotypeThreshold=0.01)
        standardGeneric("getMatrix"))
setMethod("getMatrix", signature(obj = "htestPhenStat"),
        function(obj,phenotypeThreshold=0.01){
            result <- matrix(0,2,1)
            result[1] <- sprintf("%.4f",round(pvalue(obj), digits=4))
            result[2] <-  paste(format(obj@ES, nsmall = 0),"%",sep="")
            rownames(result) <- c("p-value:","Effect size:")
            colnames(result) <- c("")
            
            if (round(pvalue(obj), digits=4) <= phenotypeThreshold){
                rownames(result)[1] <- paste("*",rownames(result)[1])
                rownames(result)[2] <- paste("*",rownames(result)[2])
            }
            else{
                rownames(result)[1] <- paste(" ",rownames(result)[1])
                rownames(result)[2] <- paste(" ",rownames(result)[2])
            }
            result
        }
)
##------------------------------------------------
setGeneric("getColumnView",
        function(obj)
        standardGeneric("getColumnView"))
setMethod("getColumnView", signature(obj = "htestPhenStat"),
        function(obj){
            return(c(subsetText(obj),sprintf("%.4f",round(pvalue(obj), digits=4)),
                            paste(format(obj@ES, nsmall = 0),"%",sep="")))
            
        }
)
##------------------------------------------------------------------------------
setGeneric("getPercentageMatrix",
        function(obj)
        standardGeneric("getPercentageMatrix"))
setMethod("getPercentageMatrix", signature(obj = "htestPhenStat"),
        function(obj){
           countMatrix <- matrixCount(obj)
           return(prop.table(countMatrix,margin=2)*100)
             
        }
 )
##------------------------------------------------------------------------------