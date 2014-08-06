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
## impress_lists.R contains functions to get and to print out IMPC objects retrieved 
## from IMPC database by using Impress SOLR REST API
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## Returns name (fieldNameTo) of the IMPC object by id (fieldValueFrom) and object class (fieldNameFrom)
getName <- function(fieldNameFrom,fieldNameTo,fieldValueFrom)
{
    fieldValueFrom <- gsub(":","\\\\:",fieldValueFrom)
    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=",fieldNameFrom,":",
                    fieldValueFrom,"&rows=0&wt=json&facet=true&"
                    ,"facet.field=",fieldNameTo,sep=""))
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    names <- unlist(json_data$facet_counts$facet_fields)
    
    numDocs <- names[seq(2,length(names),2)]
    numDocs <- as.numeric(numDocs)
    selected <- numDocs>0
    
    names <- names[seq(1,length(names),2)]
    names(names) <- NULL
    return (unlist(names[selected]))
}
##------------------------------------------------------------------------------
## Phenotyping center
getPhenCenters <- function()
{
    json_file <- "http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=*%3A*&rows=0&wt=json&facet=true&facet.field=phenotyping_center"
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    centers <- unlist(json_data$facet_counts$facet_fields$phenotyping_center)
    centers <- centers[seq(1,length(centers),2)]
    
    return (as.list(centers))
}
##------------------------------------------------------------------------------
## Phenotyping center
printPhenCenters <- function()
{
    print(unlist(getPhenCenters()))
}
##------------------------------------------------------------------------------
## Pipelines within phenotyping center
getPipelines <- function(PhenCenterName=NULL)
{
    if(is.null(PhenCenterName)){
        stop("Please define phenotyping center")
    }
    else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")
        
    }

    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
            PhenCenterName,"&rows=0&wt=json&facet=true&"
            ,"facet.field=pipeline_stable_id",sep=""))
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    pipeline_ids <- unlist(json_data$facet_counts$facet_fields$pipeline_stable_id)
    
    numDocs <- pipeline_ids[seq(2,length(pipeline_ids),2)]
    numDocs <- as.numeric(numDocs)
    #print(numDocs)
    selected <- numDocs>0
    pipeline_ids <- pipeline_ids[seq(1,length(pipeline_ids),2)]
    
    result_ids <- as.list(pipeline_ids[selected])

    return (unlist(result_ids))
   
}
##------------------------------------------------------------------------------
## Pipelines within phenotyping center
printPipelines <- function(PhenCenterName=NULL)
{
    if(is.null(PhenCenterName)){
        stop("Please define phenotyping center")
    }
    else {
        listPipelines  <- getPipelines(PhenCenterName)
        for (pipelineIndex in 1:length(listPipelines)) {
            print(paste(listPipelines[pipelineIndex],"-",
                            getName("pipeline_stable_id","pipeline_name",listPipelines[pipelineIndex])))
        }
    }
}
##------------------------------------------------------------------------------
## Procedures within pipeline of phenotyping center
getProcedures <- function(PhenCenterName=NULL, PipelineID=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)){
        stop("Please define phenotyping center and pipeline")
    }
    else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")
        
    }
    
    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                    PhenCenterName," AND pipeline_stable_id:",
                    PipelineID,"&rows=0&wt=json&facet=true&"
                    ,"facet.field=procedure_stable_id",sep=""))
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    procedures <- unlist(json_data$facet_counts$facet_fields$procedure_stable_id)
    numDocs <- procedures[seq(2,length(procedures),2)]
    numDocs <- as.numeric(numDocs)
    #print(numDocs)
    selected <- numDocs>0
    procedures <- procedures[seq(1,length(procedures),2)]
    
    return (as.list(procedures[selected]))
}
##------------------------------------------------------------------------------
## Procedures within pipeline of phenotyping center
printProcedures <- function(PhenCenterName=NULL, PipelineID=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)){
        stop("Please define phenotyping center and pipeline")
    }
    else {
        listProcedures  <- getProcedures(PhenCenterName,PipelineID)
        for (procedureIndex in 1:length(listProcedures)) {
            print(paste(listProcedures[procedureIndex],"-",
                            getName("procedure_stable_id","procedure_name",listProcedures[procedureIndex])))
        }
    }
}
##------------------------------------------------------------------------------
## Parameters measured within procedure of pipeline of phenotyping center
getParameters <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)){
        stop("Please define phenotyping center, pipeline and procedure")
    }
    else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")
        
    }
    
    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                    PhenCenterName," AND pipeline_stable_id:",
                    PipelineID," AND procedure_stable_id:",
                    ProcedureID,"&rows=0&wt=json&facet=true&"
                    ,"facet.field=parameter_stable_id",sep=""))
   # print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    parameters <- unlist(json_data$facet_counts$facet_fields$parameter_stable_id)
    numDocs <- parameters[seq(2,length(parameters),2)]
    numDocs <- as.numeric(numDocs)
    #print(numDocs)
    selected <- numDocs>0
    parameters <- parameters[seq(1,length(parameters),2)]
    
    return (as.list(parameters[selected]))
}
##------------------------------------------------------------------------------
## Parameters measured within procedure of pipeline of phenotyping center
printParameters <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)){
        stop("Please define phenotyping center, pipeline and procedure")
    }
    else {
        listParameters  <- getParameters(PhenCenterName,PipelineID,ProcedureID)
        for (parameterIndex in 1:length(listParameters)) {
            print(paste(listParameters[parameterIndex],"-",
                            getName("parameter_stable_id","parameter_name",listParameters[parameterIndex])))
        }
    }
}
##------------------------------------------------------------------------------
## Strains
getStrains <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, ParameterID=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)||is.null(ParameterID)){
        stop("Please define phenotyping center, pipeline, procedure and parameter of interest")
    }
    else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")
        
    }

    
    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                    PhenCenterName," AND pipeline_stable_id:",
                    PipelineID," AND procedure_stable_id:",
                    ProcedureID," AND parameter_stable_id:",
                    ParameterID,"&rows=0&wt=json&facet=true&"
                    ,"facet.field=strain_accession_id",sep=""))
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    strains <- unlist(json_data$facet_counts$facet_fields$strain)
    numDocs <- strains[seq(2,length(strains),2)]
    numDocs <- as.numeric(numDocs)
    #print(numDocs)
    selected <- numDocs>0
    strains <- strains[seq(1,length(strains),2)]
    
    return (as.list(strains[selected]))
}
##------------------------------------------------------------------------------
## Strains
printStrains <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, ParameterID=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)||is.null(ParameterID)){
        stop("Please define phenotyping center, pipeline, procedure and parameter of interest")
    }
    else {
        listStrains  <- getStrains(PhenCenterName,PipelineID,ProcedureID,ParameterID)
        for (strainIndex in 1:length(listStrains)) {
            print(paste(listStrains[strainIndex],"-",
                            getName("strain_accession_id","strain_name",listStrains[strainIndex])))
        }
    }
}
##------------------------------------------------------------------------------
## Genes
getGenes <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, ParameterID=NULL, StrainID=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)||is.null(ParameterID)){
        stop("Please define phenotyping center, pipeline, procedure and parameter of interest")
    }
    else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")
        
    }
    
    if (is.null(StrainID)){
    
        json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                        PhenCenterName," AND pipeline_stable_id:",
                        PipelineID," AND procedure_stable_id:",
                        ProcedureID," AND parameter_stable_id:",
                        ParameterID,"&rows=0&wt=json&facet=true&"
                        ,"facet.field=gene_accession_id",sep=""))
    }
    else {
        StrainID <- gsub(":","\\\\:",StrainID)
        json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                        PhenCenterName," AND pipeline_stable_id:",
                        PipelineID," AND procedure_stable_id:",
                        ProcedureID," AND parameter_stable_id:",
                        ParameterID," AND strain_accession_id:",
                        StrainID,"&rows=0&wt=json&facet=true&"
                        ,"facet.field=gene_accession_id",sep=""))
    }
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    genes <- unlist(json_data$facet_counts$facet_fields$gene_accession_id)
    numDocs <- genes[seq(2,length(genes),2)]
    numDocs <- as.numeric(numDocs)
    #print(numDocs)
    selected <- numDocs>0
    genes <- genes[seq(1,length(genes),2)]
    
    return (as.list(genes[selected]))
}
##------------------------------------------------------------------------------
## Genes
printGenes <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, ParameterID=NULL, StrainID=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)||is.null(ParameterID)){
        stop("Please define phenotyping center, pipeline, procedure and parameter of interest")
    }
    else {
        listGenes  <- getGenes(PhenCenterName,PipelineID,ProcedureID,ParameterID,StrainID)
        for (geneIndex in 1:length(listGenes)) {
            print(paste(listGenes[geneIndex],"-",
                            getName("gene_accession_id","gene_symbol",listGenes[geneIndex])))
        }
    }
}
##------------------------------------------------------------------------------
## Alleles
getAlleles <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, ParameterID=NULL, StrainID=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)||is.null(ParameterID)){
        stop("Please define phenotyping center, pipeline, procedure and parameter of interest")
    }
    else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")
        
    }
    add_this <- ""
    if (!is.null(StrainID)){
        StrainID <- gsub(":","\\\\:",StrainID)
        add_this <- paste(add_this," AND strain_accession_id:", StrainID, sep="")
    }  
    
    
    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                    PhenCenterName," AND pipeline_stable_id:",
                    PipelineID," AND procedure_stable_id:",
                    ProcedureID," AND parameter_stable_id:",
                    ParameterID,add_this,"&rows=0&wt=json&facet=true&"
                    ,"facet.field=allele_accession_id",sep=""))
    
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    alleles <- unlist(json_data$facet_counts$facet_fields$allele_accession_id)
    numDocs <- alleles[seq(2,length(alleles),2)]
    numDocs <- as.numeric(numDocs)
    #print(numDocs)
    selected <- numDocs>0
    alleles <- alleles[seq(1,length(alleles),2)]
    
    return (as.list(alleles[selected]))
}
##------------------------------------------------------------------------------
## Alleles
printAlleles <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, ParameterID=NULL, StrainID=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)||is.null(ParameterID)){
        stop("Please define phenotyping center, pipeline, procedure and parameter of interest")
    }
    else {
        listAlleles  <- getAlleles(PhenCenterName,PipelineID,ProcedureID,ParameterID,StrainID)
        for (alleleIndex in 1:length(listAlleles)) {
            print(paste(listAlleles[alleleIndex],"-",
                            getName("allele_accession_id","allele_symbol",listAlleles[alleleIndex])))
        }
    }
}
##------------------------------------------------------------------------------
## Zygosities
getZygosities <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, 
                          ParameterID=NULL, StrainID=NULL, GeneID=NULL, AlleleID=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)||is.null(ParameterID)){
        stop("Please define phenotyping center, pipeline, procedure and parameter of interest")
    }
    else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")
        
    }
    add_this <- ""
    if (!is.null(StrainID)){
        StrainID <- gsub(":","\\\\:",StrainID)
        add_this <- paste(add_this," AND strain_accession_id:", StrainID, sep="")
    }  
    if (!is.null(GeneID)){
        GeneID <- gsub(":","\\\\:",GeneID)
        add_this <- paste(add_this," AND gene_accession_id:", GeneID, sep="")
    } 
    if (!is.null(AlleleID)){
        AlleleID <- gsub(":","\\\\:",AlleleID)
        add_this <- paste(add_this," AND allele_accession_id:", AlleleID, sep="")
    } 
        
    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                    PhenCenterName," AND pipeline_stable_id:",
                    PipelineID," AND procedure_stable_id:",
                    ProcedureID," AND parameter_stable_id:",
                    ParameterID,add_this,"&rows=0&wt=json&facet=true&"
                    ,"facet.field=zygosity",sep=""))
    
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    zygosities <- unlist(json_data$facet_counts$facet_fields$zygosity)
    numDocs <- zygosities[seq(2,length(zygosities),2)]
    numDocs <- as.numeric(numDocs)
    #print(numDocs)
    selected <- numDocs>0
    zygosities <- zygosities[seq(1,length(zygosities),2)]
    
    return (as.list(zygosities[selected]))
}
##------------------------------------------------------------------------------
## Zygosities
printZygosities <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, 
        ParameterID=NULL, StrainID=NULL, GeneID=NULL, AlleleID=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)||is.null(ParameterID)){
        stop("Please define phenotyping center, pipeline, procedure and parameter of interest")
    }
    else {
        listZygosities  <- getZygosities(PhenCenterName,PipelineID,ProcedureID,ParameterID,StrainID,GeneID,AlleleID)
        print(unlist(listZygosities))
    }
}
##------------------------------------------------------------------------------
