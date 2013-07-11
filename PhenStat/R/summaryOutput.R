# Copyright © 2011-2013 EMBL - European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License"); 
# you may not use this file except in compliance with the License.  
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# simpleOutput.R contains simpleOutput function

summaryOutput <- function(phenTestResult,phenotypeThreshold=0.01)
# Wrapper to prepare the output of the modeling and testing results in simple user friendly form
{
    message("Model details:")
    
    message(paste("Was batch significant?",phenTestResult$model.effect.batch))
    
    message(paste("Was variance equal?",phenTestResult$model.effect.variance))
    
    if (phenTestResult$model.effect.interaction)
        sexualDimorphism = "yes"
    else 
        sexualDimorphism = "no"
    message(paste("Was there evidence of sexual dimorphism? ",sexualDimorphism," (p-value ",round(phenTestResult$model.output.interaction,digits=3),")",sep=""))
    
    message(paste("Final fitted model:",toString(phenTestResult$model.formula.genotype)))
    
    
    message("Model output:")
    
    message(paste("Genotype effect:",round(phenTestResult$model.output.genotype.nulltest.pVal,digits=3)))
    
    message(paste("Classification tag:", classificationTag(phenTestResult,phenotypeThreshold=phenotypeThreshold)))
    
    summary(phenTestResult$model.output)$tTable
    
    
}