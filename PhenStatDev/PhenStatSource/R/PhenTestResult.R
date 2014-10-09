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
## PhenTestResult.R contains PhenTestResult function
##------------------------------------------------------------------------------
PhenTestResult <- function(model.output=NULL, 
        model.dataset=NULL,
        depVariable=NULL, 
        equation="withWeight", 
        method="MM", 
        model.effect.batch=NULL, 
        model.effect.variance=NULL, 
        model.effect.sex=NULL, 
        model.effect.interaction=NULL, 
        model.output.interaction=NULL, 
        model.effect.weight=NULL, 
        numberSexes=NULL,
        pThreshold=0.05, 
        model.formula.null=NULL, 
        model.formula.genotype=NULL, 
        model.output.genotype.nulltest.pVal=NULL, 
        model.output.quality=NULL,
        model.output.summary=NULL,
        model.output.averageRefGenotype=NULL,
        model.output.percentageChanges=NULL) 

## Construct PhenTestResult object from components

{    
    x <- new("PhenTestResult",list(model.output=model.output))
    
    x$depVariable <- depVariable
    x$equation <- equation
    x$method <- method
    x$model.effect.batch <- model.effect.batch
    x$model.effect.variance <- model.effect.variance
    x$model.effect.sex <- model.effect.sex
    x$model.effect.interaction <- model.effect.interaction
    x$model.effect.weight <- model.effect.weight
    x$model.output.genotype.nulltest.pVal <- model.output.genotype.nulltest.pVal
    x$pThreshold <- pThreshold
    x$model.output <- model.output
    x$model.dataset <- model.dataset
    x$model.output.quality <- model.output.quality
    x$model.output.summary <- model.output.summary
    x$model.output.interaction <- model.output.interaction
    x$model.output.averageRefGenotype <- model.output.averageRefGenotype
    x$model.output.percentageChanges <- model.output.percentageChanges
    x
}
##------------------------------------------------------------------------------