## Copyright © 2014 EMBL - European Bioinformatics Institute
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
## impress_sets.R contains functions that are using ... REST API to retrieve datasets 
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
# R wrapper around java class getExperimentDTO from Impress REST API
getIMPCDataset <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, ParameterID=NULL, 
        AlleleID=NULL, StrainID=NULL){
    
 # Example of dataset with different metagroups:
 # df <- getIMPCDataset('WTSI','ESLIM_001','ESLIM_021_001','ESLIM_021_001_001','MGI:4362924')   
    url_main <- "http://www.mousephenotype.org/data/exportraw?"
    
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)||is.null(ParameterID)){
        stop("Please define phenotyping center, pipeline, procedure and parameter of interest")
    }
    
    add_this <- ""
    if (!is.null(StrainID)){
        add_this <- paste(add_this,"&strain_accession_id=", StrainID, sep="")
    } 
    
    url <- URLencode(paste(url_main,
                    "phenotyping_center=",PhenCenterName,
                    "&pipeline_stable_id=", PipelineID,
                    "&procedure_stable_id=", ProcedureID,
                    "&parameter_stable_id=", ParameterID,
                    "&allele_accession=", AlleleID, 
                    add_this,sep=""))
    
    #print(url)
    tryCatch(
    {
        df <- read.csv(url)
        df_trimmed <- as.data.frame(lapply(df,function (x) gsub("^\\s+|\\s+|\\n+$", "", x)))
        return (df_trimmed)
    },
    warning = function(war){
        print(paste("ERROR:  ",war))
        print(url)
    },
    error = function(err){
        print(paste("ERROR: ",err))
        print(paste("URL:",url))
    })
   
}   

# alleleAccession
# return TSV of all mutants and control discrimited

# WRAPS     public List<ExperimentDTO> getExperimentDTO(

#        // Mandatory
#        String parameterId string IMPC_001
#        String procedureId string IMPC_CBC_001
#        String pipelineId  string IMPC_CBC_002_001
#        String phenotypingCenterId string "JAX", "WTSI"
#        String alleleAccession "MGI:3452346"

#        // Optional
#        SexType sex "female"
#        List<String> zygosity "homozygous" "hemizygous" "heterozygous"
#        String strain "C57BL/6"


#        TSV RESPONSE:
#        experiment DTO unique ID
#        metadata_group
#        zygosity
#        sex
#        gene
#        allele
#        pipeline
#        procedure
#        parameter
#        phenotyping_center
#        date_of_experiment

#        // future
#        weight

#        ... other stuff

##------------------------------------------------------------------------------
## Prints or saves as csv file the table of Impress objects combinations together with the appropriate for this
## combination call of getImpressDataset function. Table can be filtered out by user. Additional table parser is needed
## in oder to call getImpressDataset multiple times (ideal case is to call function in parallel). 
getIMPCTable <- function(fileName="PhenData_IMPC",
                PhenCenterName=NULL, PipelineID=NULL, 
                ProcedureID=NULL, ParameterID=NULL, AlleleID=NULL, StrainID=NULL,
                multipleFiles=TRUE,recordsPerFile=10000)#, SexType=NULL, Zygosity=NULL){
{
        
    message(paste("Start",Sys.time()))
    if(!multipleFiles){
            message("Warning:\nAll records will be saved into one file. The size of the resulting file can be enormous.")
    }
    if(multipleFiles && recordsPerFile<1000){
            message("Warning\nMinimal number of records per file is 1000.")
            recordsPerFile <- 1000
    }
        
    if (!is.null(PhenCenterName)){
        listCenters <- c(PhenCenterName)
    } 
    else {
        listCenters <- unlist(getPhenCenters())
    }
    # ??????
    #        SexType sex "female"
    #        List<String> zygosity "homozygous" "hemizygous" "heterozygous"
    #        String strain "C57BL/6" name or id ?
    
    header <- c("Phenotyping Center","Pipeline","Screen/Procedure"
            ,"Parameter","Allele","Function to get IMPC Dataset")
    
    countRows <- 1
    totalRows <- 1
    fileNumber <- 1
    if (multipleFiles){
            currentFileName <- paste(fileName,"_",fileNumber,".csv",sep="")
    }    
    else{
            currentFileName <- fileName
    }
        
    write.table(as.matrix(t(header)), file=currentFileName, sep=",",
                col.names = FALSE,row.names = FALSE)
                
    for (centerIndex in 1:length(listCenters) ) {
        #print(listCenters[centerIndex])
        if (!is.null(PipelineID)){
            listPipelines <- c(PipelineID)
        } 
        else {
            listPipelines <- unlist(getPipelines(listCenters[centerIndex]))
        }
        pipelineIndex = 0
        while (pipelineIndex<length(listPipelines)){
            pipelineIndex <- pipelineIndex+1
        #for (pipelineIndex in 1:length(listPipelines)) {
            #print(paste(listCenters[centerIndex],listPipelines[pipelineIndex]))
            pipeline_name <- paste(
                    getName("pipeline_stable_id","pipeline_name",listPipelines[pipelineIndex]),
                    " (",listPipelines[pipelineIndex],")",sep="")
            if (!is.null(ProcedureID)){
                listProcedures <- c(ProcedureID)
            }
            else {
                listProcedures <- unlist(getProcedures(listCenters[centerIndex],listPipelines[pipelineIndex]))
            } 
            procedureIndex = 0
            while (procedureIndex<length(listProcedures)){
                procedureIndex <- procedureIndex+1
            #for (procedureIndex in 1:length(listProcedures)) {
                #print(paste(listCenters[centerIndex],listPipelines[pipelineIndex],listProcedures[procedureIndex]))
                procedure_name <- paste(
                        getName("procedure_stable_id","procedure_name",listProcedures[procedureIndex]),
                        " (",listProcedures[procedureIndex],")",sep="")
                if (!is.null(ParameterID)){
                    listParameters <- c(ParameterID)
                }
                else {
                    listParameters <- unlist(getParameters(listCenters[centerIndex],
                                    listPipelines[pipelineIndex],listProcedures[procedureIndex]))
                } 
                
                parameterIndex = 0
                while (parameterIndex<length(listParameters)){
                     parameterIndex <- parameterIndex+1
                #for (parameterIndex in 1:length(listParameters)) {
                    #print(paste(listCenters[centerIndex],listPipelines[pipelineIndex],
                    #                        listProcedures[procedureIndex],listParameters[parameterIndex]))
                    #if (!is.null(StrainID)){
                    #    listStrains <- c(StrainID)
                    #} 
                    #else {
                    #    listStrains <- unlist(getStrains(listCenters[centerIndex],
                    #                    listPipelines[pipelineIndex],listProcedures[procedureIndex],
                    #                    listParameters[parameterIndex]))
                    #}
                    #for (strainIndex in 1:length(listStrains)) {
                    #print(paste(listCenters[centerIndex],listPipelines[pipelineIndex],
                    #            listProcedures[procedureIndex],listParameters[parameterIndex],
                    #            listStrains[strainIndex]))
                    parameter_name <- paste(
                                getName("parameter_stable_id","parameter_name",listParameters[parameterIndex]),
                                " (",listParameters[parameterIndex],")",sep="")
                        
                    if (!is.null(AlleleID)){
                        listAlleles <- c(AlleleID)
                    } 
                    else {
                        listAlleles <- unlist(getAlleles(listCenters[centerIndex],
                                        listPipelines[pipelineIndex],listProcedures[procedureIndex],
                                        listParameters[parameterIndex]))#,listStrains[strainIndex]))
                    }
                    
                    alleleIndex = 0
                    while (alleleIndex<length(listAlleles)){
                        alleleIndex <- alleleIndex+1
                
                        countRows <- countRows + 1
                        totalRows <- totalRows + 1
                        #print(paste(listCenters[centerIndex],listPipelines[pipelineIndex],
                        #            listProcedures[procedureIndex],listParameters[parameterIndex],
                        #           listAlleles[alleleIndex])) 
                        function_call <- paste("getIMPCDataset('",listCenters[centerIndex],"','",
                            listPipelines[pipelineIndex],"','",
                            listProcedures[procedureIndex],"','",
                            listParameters[parameterIndex],"','",
                            listAlleles[alleleIndex],"')",sep="")
                        allele_name <- paste(
                            getName("allele_accession_id","allele_symbol",listAlleles[alleleIndex]),
                            " (",listAlleles[alleleIndex],")",sep="")
                        row_values <- c(listCenters[centerIndex],pipeline_name,
                                    procedure_name,parameter_name,
                                    allele_name, function_call)
                        if (multipleFiles && countRows>=recordsPerFile){
                            fileNumber <- fileNumber + 1
                            currentFileName <- paste(fileName,"_",fileNumber,".csv",sep="")
                            countRows <- 1
                            write.table(as.matrix(t(header)), file=currentFileName, sep=",",
                            col.names = FALSE,row.names = FALSE)
                    
                        }
                        write.table(as.matrix(t(row_values)), file = currentFileName, sep = ",", 
                        col.names = FALSE, row.names = FALSE, append=TRUE)  
                    }   
                    #}
                }
            }
        }
    }
    message(paste("End",Sys.time()))
    message(paste("Number of rows:",totalRows))
}
##------------------------------------------------------------------------------