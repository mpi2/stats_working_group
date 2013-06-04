PhenTestResult <- function(modelOutput=NULL, depVariable=NULL, equation="withWeight", batchEffect=NULL, varianceEffect=NULL, 
                  genderEffect=NULL, interactionEffect=NULL, interactionTest=NULL, weightEffect=NULL, outputLength=NULL,
                  pThreshold=0.05, model.formula=NULL, model.null=NULL, 
                  model.genotype=NULL, model=NULL, genotypeEffect=NULL, MM_fitquality=NULL) 
#    Construct PhenTestResult object from components
{
    
    
    
    x <- new("PhenTestResult",list(modelOutput=modelOutput))
    
    x$depVariable <- depVariable
    x$equation <- equation
    x$batchEffect <- batchEffect
    x$varianceEffect <- varianceEffect
    x$genderEffect <- genderEffect
    x$interactionEffect <- interactionEffect
    x$weightEffect <- weightEffect
    x$genotypeEffect <- genotypeEffect
    x$pThreshold <- pThreshold
    x$MM_fitquality <- MM_fitquality
    x$interactionTest <- interactionTest
    x$outputLength <- outputLength
    x
}
