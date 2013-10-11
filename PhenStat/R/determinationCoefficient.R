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
rsquared <- function(i)
# Coefficient of determination (R2) is a useful and intuitive tool for ascertaining 
# whether a model describes the data well: 
# it's the variance explained by the model. 
# R2 is 1 - the ratio of the unexplained variance by the model
# (residual variance, difference between observed and predicted values) over the total observed variance.
{
    #For models fit using lm
    if(class(i)=="lm") {
        Rsquared.mat=data.frame(Class=class(i),Marginal=summary(i)$r.squared,
                Conditional=NA,AIC=AIC(i)) } 
    #For models fit using lme4
    else if(inherits(i,"merMod") | class(i)=="merLmerTest") {
        #Get variance of fixed effects by multiplying coefficients by design matrix
        VarF=var(as.vector(fixef(i) %*% t(i@pp$X))) 
        #Get variance of random effects by extracting variance components
        VarRand=colSums(do.call(rbind,lapply(VarCorr(i),function(j) j[1])))
        #Get residual variance
        VarResid=attr(VarCorr(i),"sc")^2
        #Calculate marginal R-squared (fixed effects/total variance)
        Rm=VarF/(VarF+VarRand+VarResid)
        #Calculate conditional R-squared (fixed effects+random effects/total variance)
        Rc=(VarF+VarRand)/(VarF+VarRand+VarResid)
        #Bind R^2s into a matrix and return with AIC values
        Rsquared.mat=data.frame(Class=class(i),Marginal=Rm,Conditional=Rc,
                AIC=AIC(update(i,REML=F))) } 
    #For model fit using nlme  
    else if(class(i)=="lme") {
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
        #Bind R^2s into a matrix and return with AIC values
        Rsquared.mat=data.frame(Class=class(i),Marginal=Rm,Conditional=Rc,
                AIC=AIC(update(i,method="ML")))
    } else if(class(i)=="gls") {
        Rsquared = cor((i$fitted-i$residuals),i$fitted)^2
        Rsquared.mat=data.frame(Class=class(i),Marginal=Rsquared,
                Conditional=Rsquared,AIC=AIC(i)) 
    }    
    else
    { 
        stop("Function requires models of class lm, lme, mer, merMod or gls") 
    } 
    return (Rsquared.mat)
}

Cohenf <- function(phenTestResult){
# Local Effect Size
    if (phenTestResult$method=="MM") {
        R2_genotype <- rsquared(phenTestResult$model.output)$Marginal
        R2_null <- rsquared(phenTestResult$model.null)$Marginal
        Cohenf <- abs((R2_genotype-R2_null)/(1-R2_genotype))
        return(Cohenf)
    }
    else
        stop("Implemented only for the MM framework results")
}
