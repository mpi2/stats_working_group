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
#-----------------------------------------------------------------------------------
# determinationCoefficient.R contains rsquared and Cohenf functions
#-----------------------------------------------------------------------------------
rsquared.lme <- function(i) 
{
    if(class(i)=="lme"){
        
        #Get design matrix of fixed effects from model
        Fmat=model.matrix(eval(i$call$fixed)[-2],i$data)
        #Get variance of fixed effects by multiplying coefficients by design matrix
        VarF=var(as.vector(fixef(i) %*% t(Fmat)))
        #Get variance of random effects by extracting variance components
        VarRand=sum(suppressWarnings(as.numeric(nlme::VarCorr(i)[rownames(nlme::VarCorr(i))!=
                                        "Residual",1])),na.rm=T)
        #Get residual variance
        VarResid=as.numeric(nlme::VarCorr(i)[rownames(nlme::VarCorr(i))=="Residual",1])
        #Calculate marginal R-squared (fixed effects/total variance)
        Rm=VarF/(VarF+VarRand+VarResid)
        #Calculate conditional R-squared (fixed effects+random effects/total variance)
        Rc=(VarF+VarRand)/(VarF+VarRand+VarResid)
        #Bind R^2s into a matrix
        rsquared =data.frame(Class=class(i),Marginal=Rm,Conditional=Rc)
        
        return(rsquared)
    }
    else stop("Not lme")
}              
                
# The McFadden pseudo R2                
rsquared.gls <- function(phenTestResult)    
{
  return(1-(as.numeric(logLik(phenTestResult$model.output)/logLik(phenTestResult$model.null))))                  
}
                
    
Cohenf <- function(phenTestResult)
{
# Local Effect Size
# Cohen's f^2 = variance among group means/pooled within group variance
    if (phenTestResult$method=="MM") {
        if (class(phenTestResult$model.output)=="lme"){
                            R2_genotype <- rsquared.lme(phenTestResult$model.output)$Conditional
                            R2_null <- rsquared.lme(phenTestResult$model.null)$Conditional
                            Cohenf <- abs((R2_genotype-R2_null)/(1-R2_genotype))
        }
        else
            if(class(phenTestResult$model.output)=="gls"){
                            val=rsquared.gls(phenTestResult)
                            Cohenf=val/(1-val)
            }
        return(Cohenf)
    }
    else
        stop("Implemented only for the MM framework results")
}
