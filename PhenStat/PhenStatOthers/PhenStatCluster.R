PhenStatCluster<-function(phenList,i){
    variable <- as.character(phenList$dataset.stat$Variables[i])
    isContinuous <- phenList$dataset.stat$Continuous[i]
    if (!(variable %in% c("Batch","Genotype"))){
                if (isContinuous && !(variable %in% c("Weight")))
                    result <- testDataset(phenList, variable, method="MM")
                else
                if (!isContinuous){
                    result <- testDataset(phenList, variable, method="FE")
                }  
                else 
                    result <- testDataset(phenList, variable, method="MM",equation="withoutWeight")
                write(vectorOutput(result),paste("./",variable,".txt",sep=""))
    }

}    