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

# Diagnostictest.R contains testFinalModel function

testFinalModel<-function(object, result=NULL, equation=NULL, depVariable=NULL, pThreshold=0.05, keepList=NULL)
# Diagnostic test output for MM quality of fit
{
    require(nortest)
    
    # Check PhenList object
    if(is(object,"PhenList")) {
        x <- object$dataset          
    } else {
        x <- as.data.frame(object)
    }
    
    
    # Check PhenTestResult object
    if(is(result,"PhenTestResult")) {
        if (!is.null(result$model.output.quality)) return(result$model.output.quality)
        if (is.null(depVariable)) depVariable <- result$depVariable
        if (is.null(equation)) equation <- result$equation
        keep_weight <- result$model.effect.weight
        keep_gender <- result$model.effect.gender
        keep_interaction <- result$model.effect.interaction
        keep_batch <- result$model.effect.batch
        keep_equalvar <- result$model.effect.variance
        
        if (is.null(depVariable)) stop("Please define dependant variable")
        if (is.null(equation)) stop("Please define equation: 'withWeight' or 'withoutWeight'")
        if (is.null(keep_batch) || is.null(keep_equalvar) || is.null(keep_gender) || is.null(keep_interaction)) 
        stop ("Please run function 'testData' first")
        if (result$equation!=equation) stop(paste("Tests have been done with another equation: ",result$equation))
        if (result$depVariable!=depVariable) stop(paste("Tests have been done for another dependant variable: ",result$depVariable))
        
    }
    else{
        # Stop function if there are no enough needed input parameters
        if (is.null(keepList) || length(keepList)!=5) 
        stop("Please define the values for 'keepList' list, where for each effect/part of the model TRUE/FALSE value defines to keep it in the model or not: 
                'keepList=c(keepBatch,keepVariance,keepWeight,keepGender,keepInteraction)'")
        if (is.null(depVariable)) stop("Please define dependant variable")
        if (is.null(equation)) stop("Please define equation: 'withWeight' or 'withoutWeight'")
        
        keep_weight <- keepList[3]
        keep_gender <- keepList[4]
        keep_interaction <- keepList[5]
        keep_batch <- keepList[1]
        keep_equalvar <- keepList[2]
        result<-buildFinalModel(object, equation=equation, depVariable=depVariable, pThreshold=pThreshold, keepList=keepList)
    }
    
    
    a=levels(x$Genotype)
    numberofgenders=result$numberGenders
    
    if(!keep_weight && equation=="withWeight"){
        testresults=c(a[1], NA, a[2], NA, NA, NA)
        return(testresults)    
        
    }else{    
        res=resid(result$model.output)
        data_all= data.frame(x, res)
        genotype_no=length(a)
        data_all[, c("Gender", "Assay.Date")] = lapply(data_all[, c("Gender", "Assay.Date")], factor)
        No_batches=nlevels(data_all$Assay.Date)
        outputnumeric=is.numeric(modeloutput$apVar)
        
        Gp1 = subset(data_all, data_all$Genotype==a[1])
        Gp2 = subset(data_all, data_all$Genotype==a[2])
        No_Gp1=sum(is.finite(Gp1[ , depVariable]))
        No_Gp2=sum(is.finite(Gp2[ , depVariable]))
        
        if(keep_batch && No_batches >7 && outputnumeric){
            blups=ranef(modeloutput)
            blups_test= cvm.test(blups [ ,1])$p.value
            sdests = exp(attr(modeloutput$apVar, "Pars"))           #extract variance estimates
            Zbat = model.matrix(~ Assay.Date, model.frame( ~ Assay.Date, modeloutput$groups))    #create random effects design matrix
            ycov = (Zbat %*% t(Zbat)) * sdests["reStruct.Assay.Date"]^2 + diag(rep(1,nrow(modeloutput$groups))) * sdests["lSigma"]^2    #create estimated cov(y)
            Lt = chol(solve(ycov))  #Cholesky decomposition of inverse of cov(y) (see Houseman '04 eq. (2))
            rotres = Lt %*%  modeloutput$residuals[, "fixed"]    #rotated residuals
            rotated_residual_test=cvm.test(rotres)$p.value
        }else{
            blups_test=NA
            rotated_residual_test=NA
        }   
        
        
        
        if(No_Gp1>7){
            gp1_norm_res= cvm.test(Gp1$res)$p.value
        }else{
            gp1_norm_res= NA
        }    
        
        if(No_Gp2>7){
            gp2_norm_res= cvm.test(Gp2$res)$p.value
        }else{
            gp2_norm_res= NA
        }    
        
        testresults=c(a[1], gp1_norm_res, a[2], gp2_norm_res, blups_test, rotated_residual_test)
        
        
        
        return(testresults)        
        
    }
      
}
