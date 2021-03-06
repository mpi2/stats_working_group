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
## graphsDataset.R contains functions for dataset diagnostic plots:
## boxplotSexGenotype, boxplotSexGenotypeBatch, 
## scatterplotGenotypeWeigh
##------------------------------------------------------------------------------
## Raw data boxplot: split by sex and genotype
boxplotSexGenotype<-function(phenList, depVariable=NULL, 
        graphingName=NULL, outputMessages=TRUE){
    stop_message <- ""
    ## Checks
    if (is.null(depVariable)) 
    stop_message <- paste(stop_message,
            "Error:\nPlease define dependent variable 'depVariable'.\n",sep="")
    else {
        if (is.null(graphingName))
        graphingName <- depVariable
    }
    
    if(is(phenList,"PhenList")) {
        x <- getDataset(phenList)
        
        if (nchar(stop_message)==0){
            if (!(depVariable %in% colnames(x)))
            stop_message <- paste(stop_message,
                    "Error:\n",depVariable," column is missed in the dataset.",sep="")
            else {
                columnOfInterest <- x[,c(depVariable)]
                
                ## Test: depVariable is numeric 
                if(!is.numeric(columnOfInterest))
                stop_message <- paste(stop_message,
                        "Error:\n",depVariable," variable is not numeric. ",
                        "Can't create a plot based on it.",sep="")
            }       
          
        }
        
    } else {
        stop_message <- paste(stop_message,"Error:\nPlease define PhenList object.\n",sep="")
    }
          
    
    if (nchar(stop_message)>0){
        if (outputMessages){
            message(stop_message)
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt))
            stop()
        }
        else {
            stop(stop_message)
        }
    } 
    else {
        ## Plot creation
        numberofsexes <- length(levels(x$Sex))
        if (is.numeric(x[ ,depVariable]))   
          y_range <- c(min(x[ ,depVariable], na.rm=TRUE), 
                       max((x[ ,depVariable]), na.rm=TRUE))
        else
          y_range <- c(1, length(levels(x[ ,depVariable])))
        
        if(numberofsexes==2){
            Male <- subset(x, x$Sex=="Male")
            Female <- subset(x, x$Sex=="Female")      
            op <- par(mfrow=c(1,2))
            boxplot(Male[ , depVariable]~Male$Genotype, 
                    ylab=graphingName, xlab="Genotype",ylim=y_range)
            legend("topright", "Male", cex=1.3, bty="n")
            boxplot(Female[ , depVariable]~Female$Genotype, 
                    ylab=graphingName, xlab="Genotype",ylim=y_range)
            legend("topright", "Female", cex=1.3, bty="n")
            par(op) 
            op_normal <- par(mfrow=c(1,1))
            par(op_normal) 
        }else{
            op <- par(mfrow=c(1,1))
            boxplot(x[ ,depVariable]~x$Genotype, ylab=graphingName, xlab="Genotype")  
            par(op)  
        }
    }
}    
##------------------------------------------------------------------------------
##Graphing after accounting for batch
boxplotSexGenotypeBatchAdjusted<-function(phenList, depVariable=NULL, 
        graphingName=NULL, outputMessages=TRUE){
    stop_message <- ""
    ## Checks
    if (is.null(depVariable)) 
    stop_message <- paste(stop_message,
            "Error:\nPlease define dependent variable 'depVariable'.\n",sep="")
    else {
        if (is.null(graphingName))
        graphingName <- paste(depVariable," adjusted for batch",sep="")
    }
    
    if(is(phenList,"PhenList")) {
        x <- getDataset(phenList)
        
        if (nchar(stop_message)==0){
            if (!(depVariable %in% colnames(x)))
                stop_message <- paste(stop_message,
                    "Error:\n",depVariable," column is missed in the dataset.",sep="")
            else {
                if (!multipleBatches(phenList)){ 
                    stop_message <- paste(stop_message,
                            "Error:\n Plot creation relies on having batch variation.",sep="")
                }
                else { 
                    columnOfInterest <- x[,c(depVariable)]
                    
                    ## Test: depVariable is numeric 
                    if(!is.numeric(columnOfInterest))
                    stop_message <- paste(stop_message,
                            "Error:\n",depVariable," variable is not numeric. ",
                            "Can't create a plot based on it.",sep="")
                }
            }       
            
        }
        
    } else {
        stop_message <- paste(stop_message,"Error:\nPlease define PhenList object.\n",sep="")
    }
    
    
    if (nchar(stop_message)>0){
        if (outputMessages){
            message(stop_message)
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt))
            stop()
        }
        else {
            stop(stop_message)
        }
    } 
    else {
        message(paste("Information:\n",depVariable," variable is adjusted ",
                        "treating batch as a random effect.",sep=""))
        ## Plot creation             
        # relies on having batch variation   
        x[, depVariable] <- getColumnBatchAdjusted(phenList,depVariable)
        
        ## Plot creation
        numberofsexes <- length(levels(x$Sex))
        if (is.numeric(x[ ,depVariable]))   
          y_range <- c(min(x[ ,depVariable], na.rm=TRUE), 
                       max((x[ ,depVariable]), na.rm=TRUE))
        else
          y_range <- c(1, length(levels(x[ ,depVariable])))
        
        if(numberofsexes==2){
            Male <- subset(x, x$Sex=="Male")
            Female <- subset(x, x$Sex=="Female")            
            op <- par(mfrow=c(1,2))
            boxplot(Male[ , depVariable]~Male$Genotype, 
                    ylab=graphingName, xlab="Genotype", ylim=y_range)
            legend("topright", "Male", cex=1.3, bty="n")
            boxplot(Female[ , depVariable]~Female$Genotype, 
                    ylab=graphingName, xlab="Genotype", ylim=y_range)
            legend("topright", "Female", cex=1.3, bty="n")
            par(op) 
            op_normal <- par(mfrow=c(1,1))
            par(op_normal) 
        }else{
            op <- par(mfrow=c(1,1))
            boxplot(x[ ,depVariable]~x$Genotype, ylab=graphingName, xlab="Genotype")  
            par(op)  
        }
        

    }
}
##------------------------------------------------------------------------------
##Graphing after accounting for weight and batch
boxplotSexGenotypeWeightBatchAdjusted<-function(phenList, depVariable=NULL, 
                                          graphingName=NULL, outputMessages=TRUE){
  stop_message <- ""
  ## Checks
  if (is.null(depVariable)) 
    stop_message <- paste(stop_message,
                          "Error:\nPlease define dependent variable 'depVariable'.\n",sep="")
  else {
    if (is.null(graphingName))
      graphingName <- paste(depVariable," adjusted for weight and batch",sep="")
  }
  
  if(is(phenList,"PhenList")) {
    x <- getDataset(phenList)
    
    if (nchar(stop_message)==0){
      if (!(depVariable %in% colnames(x)))
        stop_message <- paste(stop_message,
          "Error:\n",depVariable," column is missed in the dataset.",
          sep="")
      else {
        if (!multipleBatches(phenList) && !weightIn(phenList)){ 
          stop_message <- paste(stop_message,
          "Error:\n Plot creation relies on having batch and weight variation.",
          sep="")
        }
        else { 
          columnOfInterest <- x[,c(depVariable)]
          
          ## Test: depVariable is numeric 
          if(!is.numeric(columnOfInterest))
            stop_message <- paste(stop_message,
                            "Error:\n",depVariable," variable is not numeric. ",
                            "Can't create a plot based on it.",sep="")
        }
      }       
      
    }
    
  } else {
    stop_message <- paste(stop_message,
                          "Error:\nPlease define PhenList object.\n",sep="")
  }
  
  
  if (nchar(stop_message)>0){
    if (outputMessages){
      message(stop_message)
      opt <- options(show.error.messages=FALSE)
      on.exit(options(opt))
      stop()
    }
    else {
      stop(stop_message)
    }
  } 
  else {
    message(paste("Information:\n",depVariable," variable is adjusted ",
            "treating batch as a random effect and weight as a fixed effect.",
            sep=""))
    ## Plot creation             
    # relies on having batch variation   
    x[, depVariable] <- getColumnWeightBatchAdjusted(phenList,depVariable)
    
    ## Plot creation
    numberofsexes <- length(levels(x$Sex))
    if (is.numeric(x[ ,depVariable]))   
      y_range <- c(min(x[ ,depVariable], na.rm=TRUE), 
                   max((x[ ,depVariable]), na.rm=TRUE))
    else
      y_range <- c(1, length(levels(x[ ,depVariable])))
    
    if(numberofsexes==2){
      Male <- subset(x, x$Sex=="Male")
      Female <- subset(x, x$Sex=="Female")
      op <- par(mfrow=c(1,2))
      boxplot(Male[ , depVariable]~Male$Genotype, 
              ylab=graphingName, xlab="Genotype", ylim=y_range)
      legend("topright", "Male", cex=1.3, bty="n")
      boxplot(Female[ , depVariable]~Female$Genotype, 
              ylab=graphingName, xlab="Genotype",ylim=y_range)
      legend("topright", "Female", cex=1.3, bty="n")
      par(op) 
      op_normal <- par(mfrow=c(1,1))
      par(op_normal) 
    }else{
      op <- par(mfrow=c(1,1))
      boxplot(x[ ,depVariable]~x$Genotype, ylab=graphingName, xlab="Genotype")  
      par(op)  
    }
    
    
  }
}
##------------------------------------------------------------------------------
## Raw data boxplot: split by sex,genotype and batch
## NB! Deprecated function  
boxplotSexGenotypeBatch<-function(phenList, depVariable=NULL, 
        graphingName=NULL, outputMessages=TRUE){
    
    message("NB! The function is deprecated. Please use scatterplotSexGenotypeBatch instead.")
    stop_message <- ""
    ## Checks
    
    if (is.null(depVariable)) 
        stop_message <- paste(stop_message,
            "Error:\nPlease define dependent variable 'depVariable'.\n",sep="")
    else {
        if (is.null(graphingName))
            graphingName <- depVariable
    }
    
    if(is(phenList,"PhenList")) {
        x <- getDataset(phenList)
        refGenotype <- refGenotype(phenList)   
        
        if (nchar(stop_message)==0){
            if (!(depVariable %in% colnames(x)))
            stop_message <- paste(stop_message,
                    "Error:\n",depVariable," column is missed in the dataset.",sep="")
            else {
                columnOfInterest <- x[,c(depVariable)]
                
                ## Test: depVariable is numeric 
                if(!is.numeric(columnOfInterest))
                stop_message <- paste(stop_message,
                        "Error:\n",depVariable," variable is not numeric. ",
                        "Can't create a plot based on it.",sep="")
            }       
            
            
            if (!('Batch' %in% colnames(x)))
            stop_message <- paste(stop_message,
                    "Error:\nBatch column is missed in the dataset.\n",sep="")
        }
             
    } else {
        stop_message <- paste(stop_message,"Error:\nPlease define PhenList object.\n",sep="")
    }

            
   
    
    
    
    
    if (nchar(stop_message)>0){
        if (outputMessages){
            message(stop_message)
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt))
            stop()
        }
        else {
            stop(stop_message)
        }
    } 
    else {
        ## Plot creation
        numberofsexes <- length(levels(x$Sex))
        if (is.numeric(x[ ,depVariable]))   
          y_range <- c(min(x[ ,depVariable], na.rm=TRUE), 
                  max((x[ ,depVariable]), na.rm=TRUE))
        else
          y_range <- c(1, length(levels(x[ ,depVariable])))
        
        if(numberofsexes==2){
            Male <- subset(x, x$Sex=="Male")
            Female <- subset(x, x$Sex=="Female")
            Male$Batch <- factor(Male$Batch)
            Female$Batch <- factor(Female$Batch)
            
            op <- par(mfrow=c(1,2)) 
            
            boxplot(Male[ , depVariable]~Male$Genotype+Male$Batch, 
                    subset=(Male$Genotype==refGenotype), 
                    ylab=graphingName, ylim=y_range, xlab="Batch",  names=NULL)
            
            boxplot(Male[ , depVariable]~Male$Genotype + Male$Batch, add=TRUE, 
                    subset=(Male$Genotype!=refGenotype), ylim=y_range, 
                    ylab=graphingName, xlab="Batch",  col="red", names=NULL)
            
            legend("topright", "Male", cex=1.3, bty="n")
            
            boxplot(Female[ , depVariable]~Female$Genotype + Female$Batch, 
                    subset=(Female$Genotype==refGenotype),ylim=y_range,
                    ylab=graphingName, xlab="Batch",  names=NULL )
            
            boxplot(Female[ , depVariable]~Female$Genotype + Female$Batch, add=TRUE, 
                    subset=(Female$Genotype!=refGenotype),ylim=y_range, 
                    ylab=graphingName, xlab="Batch",  col="red",names=NULL)
            
            legend("topright", "Female", cex=1.3, bty="n")
            
            par(op)
            
            op_normal <- par(mfrow=c(1,1))
            par(op_normal) 
        }else{
            op <- par(mfrow=c(1,1))
            boxplot(x[ ,depVariable]~x$Genotype+x$Batch,
                    subset=(x$Genotype==refGenotype), 
                    ylab=graphingName, xlab="Batch", names=NULL) 
            
            boxplot(x[ ,depVariable]~x$Genotype+x$Batch,
                    subset=(x$Genotype!=refGenotype), add=TRUE, 
                    ylab=graphingName, xlab="Batch",  col="red", names=NULL)    
            
            par(op)
        }
    }    
}    
##------------------------------------------------------------------------------
## Raw data boxplot: split by sex,genotype and batch 
scatterplotSexGenotypeBatch<-function(phenList, depVariable=NULL,        
        graphingName=NULL, outputMessages=TRUE){    
    stop_message <- ""
    
    ## Checks    
    if (is.null(depVariable))
        stop_message <- paste(stop_message,
            "Error:\nPlease define dependent variable 'depVariable'.\n",sep="")
    else {
        if (is.null(graphingName))
        graphingName <- depVariable
    }
    if(is(phenList,"PhenList")) {
        x <- getDataset(phenList)
        refGenotype <- refGenotype(phenList)
        if (nchar(stop_message)==0){
            if (!(depVariable %in% colnames(x)))
                stop_message <- paste(stop_message,
                    "Error:\n",depVariable," column is missed in the dataset.",sep="")
            else {
                columnOfInterest <- x[,c(depVariable)]
                ## Test: depVariable is numeric
                if(!is.numeric(columnOfInterest))
                    stop_message <- paste(stop_message,
                        "Error:\n",depVariable," variable is not numeric. ",
                        "Can't create a plot based on it.",sep="")
            }      
            if (!('Batch' %in% colnames(x)))
                stop_message <- paste(stop_message,
                    "Error:\nBatch column is missed in the dataset.\n",sep="")
        }
        
    } else {
        stop_message <- paste(stop_message,"Error:\nPlease define PhenList object.\n",sep="")
    }
    
    if (nchar(stop_message)>0){
        if (outputMessages){
            message(stop_message)
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt))
            stop()
        }
        else {
            stop(stop_message)
        }
    }
    else {
        ## Plot creation
        numberofsexes <- length(levels(x$Sex))
        if (is.numeric(x[ ,depVariable]))  
            y_range <- c(min(x[ ,depVariable], na.rm=TRUE),
                max((x[ ,depVariable]), na.rm=TRUE))
        else
            y_range <- c(1, length(levels(x[ ,depVariable])))
        
       if(numberofsexes==2){
            Male <- subset(x, x$Sex=="Male")
            Female <- subset(x, x$Sex=="Female")
            Male$Batch <- factor(Male$Batch)
            Female$Batch <- factor(Female$Batch)
            op <- par(mfrow=c(1,2))
            stripchart(Male[ , depVariable]~Male$Batch, ,pch=1,vertical=T, 
                    subset=(Male$Genotype==refGenotype),
                    ylab=graphingName, ylim=y_range, xlab="Batch", xaxt='n')
            points(Male[ , depVariable]~Male$Batch, subset=(Male$Genotype!=refGenotype), col="red")
            legend("topright", "Male", cex=1.3, bty="n")
            stripchart(Female[ , depVariable]~Female$Batch, ,pch=1,vertical=T, 
                    subset=(Female$Genotype==refGenotype), 
                    ylab=graphingName, ylim=y_range, xlab="Batch", xaxt='n')
            points(Female[ , depVariable]~Female$Batch, subset=(Female$Genotype!=refGenotype), col="red")
            legend("topright", "Female", cex=1.3, bty="n")
            par(op)
            op_normal <- par(mfrow=c(1,1))
            par(op_normal)
        }else{
            op <- par(mfrow=c(1,1))
            stripchart(x[ , depVariable]~x$Batch, ,pch=1,vertical=T, subset=(x$Genotype==refGenotype),
                    ylab=graphingName, ylim=y_range, xlab="Batch", xaxt='n')
            points(x[ , depVariable]~x$Batch, subset=(x$Genotype!=refGenotype), col="red")
            par(op)
        }
    }   
}  
##------------------------------------------------------------------------------
## Raw data scatterplot: body weight versus dependant variable
scatterplotGenotypeWeight<-function(phenList, depVariable=NULL, 
        graphingName=NULL, outputMessages=TRUE){
    stop_message <- ""
    
    
    if (is.null(depVariable)) 
    stop_message <- paste(stop_message,
            "Error:\nPlease define dependent variable 'depVariable'.\n",sep="")
    else {
        if (is.null(graphingName))
        graphingName <- depVariable
    }
    
    ## Checks
    if(is(phenList,"PhenList")) {
        x <- getDataset(phenList)     
        
        if (nchar(stop_message)==0){ 
            if (!(depVariable %in% colnames(x)))
            stop_message <- paste(stop_message,
                    "Error:\n",depVariable," column is missed in the dataset.",sep="")
            else {
                columnOfInterest <- x[,c(depVariable)]        
                ## Test: depVariable is numeric 
                if(!is.numeric(columnOfInterest))
                stop_message <- paste(stop_message,
                        "Error:\n",depVariable," variable is not numeric. ",
                        "Can't create a plot based on it.",sep="")
                
                # Checks of Weight
                checkWeight <- columnChecks(x,"Weight",4) 
                if (! checkWeight[1])
                stop_message <- paste(stop_message,
                        "Error:\nWeight column is not present in dataset. ",
                        "Can't create a plot.\n", sep="") 
                
                else{
                    if (! checkWeight[2])
                    stop_message <- paste(stop_message,
                            "Error:\nWeight column is not numeric. ",
                            "Can't create a plot based on it.\n",sep="")
                    
                    if (! checkWeight[3])
                    stop_message <- paste(stop_message,
                            "Error:\nWeight column does not have enough ",
                            "data points for for genotype/sex combinations.\n",sep="")
                    
                }
            }
        }   
        
    } else {
        stop_message <- paste(stop_message,"Error:\nPlease define PhenList object.\n",sep="")
    }

       
    

    if (nchar(stop_message)>0){
        if (outputMessages){
            message(stop_message)
            opt <- options(show.error.messages=FALSE)
            on.exit(options(opt))
            stop()
        }
        else {
            stop(stop_message)
        }
    } 
    else {
        ## Plot creation
        model.formula <- as.formula(paste(depVariable, "~", 
                        paste("Weight", "Genotype", sep= "|")))
        
        scatterplot(data=x, model.formula, ylab=graphingName, 
                pch=c(1,1),col=c("black","red"))
    }
}
##------------------------------------------------------------------------------