## Copyright © 2011-2014 EMBL - European Bioinformatics Institute
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
## lists.R contains functions that are using Impress SOLR REST API to retrieve data 
##------------------------------------------------------------------------------
library("rjson")
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
    pipelines <- unlist(json_data$facet_counts$facet_fields$pipeline_stable_id)
    numDocs <- pipelines[seq(2,length(pipelines),2)]
    numDocs < -as.numeric(numDocs)
    #print(numDocs)
    selected <- numDocs>0
    pipelines <- pipelines[seq(1,length(pipelines),2)]
    
    return (as.list(pipelines[selected]))
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
    numDocs < -as.numeric(numDocs)
    #print(numDocs)
    selected <- numDocs>0
    procedures <- procedures[seq(1,length(procedures),2)]
    
    return (as.list(procedures[selected]))
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
    numDocs < -as.numeric(numDocs)
    #print(numDocs)
    selected <- numDocs>0
    parameters <- parameters[seq(1,length(parameters),2)]
    
    return (as.list(parameters[selected]))
}
##------------------------------------------------------------------------------
## Strains used during measurments of patameter in phenotyping center
getStrains <- function(PhenCenterName=NULL, ParameterID=NULL)
{
    if(is.null(PhenCenterName)||is.null(ParameterID)){
        stop("Please define phenotyping center and parameter of interest")
    }
    else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")
        
    }
    
    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                    PhenCenterName," AND parameter_stable_id:",
                    ParameterID,"&rows=0&wt=json&facet=true&"
                    ,"facet.field=strain",sep=""))
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    strains <- unlist(json_data$facet_counts$facet_fields$strain)
    numDocs <- strains[seq(2,length(strains),2)]
    numDocs < -as.numeric(numDocs)
    #print(numDocs)
    selected <- numDocs>0
    strains <- strains[seq(1,length(strains),2)]
    
    return (as.list(strains[selected]))
}
##------------------------------------------------------------------------------
## Genes used during measurments of patameter in phenotyping center
getGenes <- function(PhenCenterName=NULL, ParameterID=NULL, strain=NULL)
{
    if(is.null(PhenCenterName)||is.null(ParameterID)){
        stop("Please define phenotyping center and parameter of interest")
    }
    else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")
        
    }
    
    if (is.null(strain)){
    
    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                    PhenCenterName," AND parameter_stable_id:",
                    ParameterID,"&rows=0&wt=json&facet=true&"
                    ,"facet.field=gene_accession",sep=""))
    }
    else {
        strain <- paste("\"",strain,"\"",sep="")
        json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                        PhenCenterName," AND parameter_stable_id:",
                        ParameterID," AND strain:",
                        strain,"&rows=0&wt=json&facet=true&"
                        ,"facet.field=gene_accession",sep=""))
    }
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    genes <- unlist(json_data$facet_counts$facet_fields$gene_accession)
    numDocs <- genes[seq(2,length(genes),2)]
    numDocs < -as.numeric(numDocs)
    #print(numDocs)
    selected <- numDocs>0
    genes <- genes[seq(1,length(genes),2)]
    
    return (as.list(genes[selected]))
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

