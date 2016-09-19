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
library("rjson")

##------------------------------------------------------------------------------
## unwrapSolrPivotResults - Function to unwrap the facet results from a solr call
unwrapSolrPivotResults <- function(facets)
{

	if (length(facets)==0) {
		return (list())
	}

	# facets is an array that looks like
	# [1] "hemizygote"   "0"  "heterozygote" "0"  "homozygote"   "0"
	# numDocs is every other value in the array
	numDocs <- facets[seq(2,length(facets),2)]
	numDocs <- as.numeric(numDocs)

	# Return only the values that have results
	selected <- numDocs>0
	results <- facets[seq(1,length(facets),2)]
	#names(results) <- NULL
	return (as.list(results[selected]))

}


##------------------------------------------------------------------------------
## Returns name (fieldNameTo) of the IMPC object by id (fieldValueFrom) and object class (fieldNameFrom)
getName <- function(fieldNameFrom,fieldNameTo,fieldValueFrom)
{
    fieldValueFrom <- gsub(":","\\\\:",fieldValueFrom)
    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=",fieldNameFrom,":",
                    fieldValueFrom,"&rows=0&wt=json&facet=true&facet.mincount=1&facet.limit=-1&"
                    ,"facet.field=",fieldNameTo,sep=""))
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    names <- unlist(json_data$facet_counts$facet_fields)
	return (unwrapSolrPivotResults(names))
}

##------------------------------------------------------------------------------
## Phenotyping center
getPhenCenters <- function(excludeLegacyPipelines=TRUE)
{
    json_file <- paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=*%3A*&rows=0&",
            "wt=json&facet=true&facet.mincount=1&facet.limit=-1&facet.field=phenotyping_center",sep="")
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    centers <- unlist(json_data$facet_counts$facet_fields$phenotyping_center)

	if (length(centers)==0) {
		return (list())
	}

    centers <- centers[seq(1,length(centers),2)]

    if (excludeLegacyPipelines){
        for (centerIndex in 1:length(centers) ) {
            listPipelines <- getPipelines(centers[centerIndex],excludeLegacyPipelines)
            if (length(listPipelines)==0){
                centers <- centers[-centerIndex]
            }
        }
    }
    return (as.list(centers))
}
##------------------------------------------------------------------------------
## Phenotyping center
printPhenCenters <- function(n=NULL, excludeLegacyPipelines=TRUE)
{
    if (is.null(n) || n>length(getPhenCenters(excludeLegacyPipelines))){
        n <- length(getPhenCenters(excludeLegacyPipelines))
    }
    print(unlist(getPhenCenters(excludeLegacyPipelines))[c(1:n)])
}
##------------------------------------------------------------------------------
## Pipelines within phenotyping center
getPipelines <- function(PhenCenterName=NULL,excludeLegacyPipelines=TRUE)
{
    if(is.null(PhenCenterName)){
        stop("Please define phenotyping center")
    } else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")

    }

    legacy_pipelines <- c("M-G-P_001","ESLIM_001","ESLIM_002","ESLIM_003","GMC_001")

    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
            PhenCenterName,"&rows=0&wt=json&facet=true&facet.mincount=1&facet.limit=-1&"
            ,"facet.field=pipeline_stable_id",sep=""))
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    pipeline_ids <- unlist(json_data$facet_counts$facet_fields$pipeline_stable_id)

	result_ids <- unwrapSolrPivotResults(pipeline_ids)

    if (excludeLegacyPipelines){
        result_ids <- result_ids[!(result_ids %in% legacy_pipelines)]
    }

    return (unlist(result_ids))

}
##------------------------------------------------------------------------------
## Pipelines within phenotyping center
printPipelines <- function(PhenCenterName=NULL, n=NULL, excludeLegacyPipelines=TRUE)
{
    if(is.null(PhenCenterName)){
        stop("Please define phenotyping center")
    } else {
        listPipelines  <- getPipelines(PhenCenterName,excludeLegacyPipelines)
        if (is.null(n) || n>length(listPipelines)){
            n <- length(listPipelines)
        }

        for (pipelineIndex in 1:n) {
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
    } else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")

    }

    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                    PhenCenterName," AND pipeline_stable_id:",
                    PipelineID,"&rows=0&wt=json&facet=true&facet.mincount=1&facet.limit=-1&"
                    ,"facet.field=procedure_stable_id",sep=""))
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    procedures <- unlist(json_data$facet_counts$facet_fields$procedure_stable_id)

	return (unwrapSolrPivotResults(procedures))
}
##------------------------------------------------------------------------------
## Procedures within pipeline of phenotyping center
printProcedures <- function(PhenCenterName=NULL, PipelineID=NULL, n=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)){
        stop("Please define phenotyping center and pipeline")
    } else {
        listProcedures  <- getProcedures(PhenCenterName,PipelineID)
        if (is.null(n) || n>length(listProcedures)){
            n <- length(listProcedures)
        }
        for (procedureIndex in 1:n) {
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
    } else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")

    }

    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                    PhenCenterName," AND pipeline_stable_id:",
                    PipelineID," AND procedure_stable_id:",
                    ProcedureID,"&rows=0&wt=json&facet=true&facet.mincount=1&facet.limit=-1&"
                    ,"facet.field=parameter_stable_id",sep=""))
   # print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    parameters <- unlist(json_data$facet_counts$facet_fields$parameter_stable_id)

	return (unwrapSolrPivotResults(parameters))
}
##------------------------------------------------------------------------------
## Parameters measured within procedure of pipeline of phenotyping center
printParameters <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, n=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)){
        stop("Please define phenotyping center, pipeline and procedure")
    } else {
        listParameters  <- getParameters(PhenCenterName,PipelineID,ProcedureID)
        if (is.null(n) || n>length(listParameters)){
            n <- length(listParameters)
        }
        for (parameterIndex in 1:n) {
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
    } else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")

    }


    json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                    PhenCenterName," AND pipeline_stable_id:",
                    PipelineID," AND procedure_stable_id:",
                    ProcedureID," AND parameter_stable_id:",
                    ParameterID,"&rows=0&wt=json&facet=true&facet.mincount=1&facet.limit=-1&"
                    ,"facet.field=strain_accession_id",sep=""))
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    strains <- unlist(json_data$facet_counts$facet_fields$strain)

	return (unwrapSolrPivotResults(strains))

}
##------------------------------------------------------------------------------
## Strains
printStrains <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, ParameterID=NULL, n=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)||is.null(ParameterID)){
        stop("Please define phenotyping center, pipeline, procedure and parameter of interest")
    } else {
        listStrains  <- getStrains(PhenCenterName,PipelineID,ProcedureID,ParameterID)
        if (is.null(n) || n>length(listStrains)){
            n <- length(listStrains)
        }
        for (strainIndex in 1:n) {
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
    } else {
        PhenCenterName <- paste("\"",PhenCenterName,"\"",sep="")

    }

    if (is.null(StrainID)){

        json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                        PhenCenterName," AND pipeline_stable_id:",
                        PipelineID," AND procedure_stable_id:",
                        ProcedureID," AND parameter_stable_id:",
                        ParameterID,"&rows=0&wt=json&facet=true&facet.mincount=1&facet.limit=-1&"
                        ,"facet.field=gene_accession_id",sep=""))
    } else {
        StrainID <- gsub(":","\\\\:",StrainID)
        json_file <- URLencode(paste("http://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=phenotyping_center:",
                        PhenCenterName," AND pipeline_stable_id:",
                        PipelineID," AND procedure_stable_id:",
                        ProcedureID," AND parameter_stable_id:",
                        ParameterID," AND strain_accession_id:",
                        StrainID,"&rows=0&wt=json&facet=true&facet.mincount=1&facet.limit=-1&"
                        ,"facet.field=gene_accession_id",sep=""))
    }
    #print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    genes <- unlist(json_data$facet_counts$facet_fields$gene_accession_id)

	return (unwrapSolrPivotResults(genes))

}
##------------------------------------------------------------------------------
## Genes
printGenes <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, ParameterID=NULL, StrainID=NULL,
        n=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)||is.null(ParameterID)){
        stop("Please define phenotyping center, pipeline, procedure and parameter of interest")
    } else {
        listGenes  <- getGenes(PhenCenterName,PipelineID,ProcedureID,ParameterID,StrainID)
        if (is.null(n) || n>length(listGenes)){
            n <- length(listGenes)
        }
        for (geneIndex in 1:n) {
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
    } else {
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
                    ParameterID,add_this,"&rows=0&wt=json&facet=true&facet.mincount=1&facet.limit=-1&"
                    ,"facet.field=allele_accession_id",sep=""))

    print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    alleles <- unlist(json_data$facet_counts$facet_fields$allele_accession_id)

	return (unwrapSolrPivotResults(alleles))

}
##------------------------------------------------------------------------------
## Alleles
printAlleles <- function(PhenCenterName=NULL, PipelineID=NULL, ProcedureID=NULL, ParameterID=NULL, StrainID=NULL,
        n=NULL)
{
    if(is.null(PhenCenterName)||is.null(PipelineID)||is.null(ProcedureID)||is.null(ParameterID)){
        stop("Please define phenotyping center, pipeline, procedure and parameter of interest")
    } else {
        listAlleles  <- getAlleles(PhenCenterName,PipelineID,ProcedureID,ParameterID,StrainID)
        if (is.null(n) || n>length(listAlleles)){
            n <- length(listAlleles)
        }
        for (alleleIndex in 1:n) {
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
    } else {
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
                    ParameterID,add_this,"&rows=0&wt=json&facet=true&facet.mincount=1&facet.limit=-1&"
                    ,"facet.field=zygosity",sep=""))

    print(json_file)
    json_data <- fromJSON(paste(readLines(json_file), collapse=""))
    zygosities <- unlist(json_data$facet_counts$facet_fields$zygosity)

	return (unwrapSolrPivotResults(zygosities))
}
##------------------------------------------------------------------------------
