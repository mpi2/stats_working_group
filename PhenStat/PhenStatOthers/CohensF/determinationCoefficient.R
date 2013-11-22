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
rsquared <- function(i,base_level_model=NULL) 
{
    rsquared = NULL
    
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
        
       
    }
    else if(class(i)=="gls" && !(is.null(base_level_model))){
        # The McFadden pseudo R2 
        rsquared = 1-(as.numeric(logLik(i)/logLik(base_level_model)))   
    }   
    
    return(rsquared)
}              
#-----------------------------------------------------------------------------------                
Cohenf.Conditional <- function(phenTestResult)
{
# Local Effect Size
# Cohen's f^2 = variance among group means/pooled within group variance
    if (phenTestResult$method=="MM") {
        if (class(phenTestResult$model.output)=="lme"){
                            R2_genotype <- rsquared(phenTestResult$model.output)$Conditional
                            R2_null <- rsquared(phenTestResult$model.null)$Conditional
                            Cohenf <- abs((R2_genotype-R2_null)/(1-R2_genotype))
        }
        else
            if(class(phenTestResult$model.output)=="gls"){
                intercept_only_model = do.call("gls", args = list(as.formula(paste(phenTestResult$depVariable," ~ 1")),  phenTestResult$model.dataset, na.action="na.exclude"))    
                R2_genotype <- rsquared(phenTestResult$model.output,intercept_only_model)
                R2_null <- rsquared(phenTestResult$model.null,intercept_only_model)
                Cohenf <- abs((R2_genotype-R2_null)/(1-R2_genotype))
            }
        return(Cohenf)
    }
}
#-----------------------------------------------------------------------------------
Cohenf.Marginal<- function(phenTestResult)
{
    # Local Effect Size
    # Cohen's f^2 = variance among group means/pooled within group variance
    if (phenTestResult$method=="MM") {
        if (class(phenTestResult$model.output)=="lme"){
            R2_genotype <- rsquared(phenTestResult$model.output)$Marginal
            R2_null <- rsquared(phenTestResult$model.null)$Marginal
            Cohenf <- abs((R2_genotype-R2_null)/(1-R2_genotype))
        }
        else
        if(class(phenTestResult$model.output)=="gls"){
            intercept_only_model = do.call("gls", args = list(as.formula(paste(phenTestResult$depVariable," ~ 1")),  phenTestResult$model.dataset, na.action="na.exclude"))    
            R2_genotype <- rsquared(phenTestResult$model.output,intercept_only_model)
            R2_null <- rsquared(phenTestResult$model.null,intercept_only_model)
            Cohenf <- abs((R2_genotype-R2_null)/(1-R2_genotype))
        }
        return(Cohenf)
    }
}
#-----------------------------------------------------------------------------------