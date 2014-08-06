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
## Diagnostictest.R contains testFinalModel function
##------------------------------------------------------------------------------
## Diagnostic test output for MM quality of fit. There are no arguments checks
## assuming that function is called internally from the finalModel function.
## Otherwise should be used with precaution.
testFinalModel<-function(phenTestResult)
{
    result <- phenTestResult
    x <- result$model.dataset
    depVariable <- result$depVariable
    equation <- result$equation
    keep_weight <- result$model.effect.weight
    keep_sex <- result$model.effect.sex
    keep_interaction <- result$model.effect.interaction
    keep_batch <- result$model.effect.batch
    keep_equalvar <- result$model.effect.variance
    
    a <- levels(x$Genotype)
    numberofsexes <- result$numberSexes
    
    if(!keep_weight && equation=="withWeight"){
        testresults <- c(a[1], NA, a[2], NA, NA, NA)
        
    }else if ('Batch' %in% colnames(x)){
        res <- resid(result$model.output)
        data_all <- data.frame(x, res)
        genotype_no <- length(a)
        data_all[, c("Sex", "Batch")] <- lapply(data_all[, c("Sex", 
                                "Batch")], factor)
        No_batches <- nlevels(data_all$Batch)
        outputnumeric <- is.numeric(result$model.output$apVar)
        Gp1 <- subset(data_all, data_all$Genotype==a[1])
        Gp2 <- subset(data_all, data_all$Genotype==a[2])
        No_Gp1 <- sum(is.finite(Gp1[ , c("res")]))
        No_Gp2 <- sum(is.finite(Gp2[ , c("res")]))
        
        if(keep_batch && No_batches >7 && outputnumeric){
            
            blups <- ranef(result$model.output)
            blups_test <- suppressWarnings(cvm.test(blups [ ,1])$p.value)
            sdests <- exp(attr(result$model.output$apVar, 
                            "Pars"))           #extract variance estimates
            
            Zbat <- model.matrix(~ Batch, model.frame( ~ Batch, 
                            result$model.output$groups))    #create random effects design matrix
            
            ycov <- (Zbat %*% t(Zbat)) * sdests["reStruct.Batch"]^2 + diag(rep
                    (1,nrow(result$model.output$groups))) * sdests["lSigma"]^2  #create estimated cov(y)  
            
            Lt <- chol(solve(ycov)) #Cholesky decomposition of inverse of cov
            # (y) (see Houseman '04 eq. (2))
            rotres <- Lt %*%  result$model.output$residuals[,"fixed"] #rotated residuals
            
            rotated_residual_test <- suppressWarnings(cvm.test(rotres)$p.value)
        }else{
            blups_test <- NA
            rotated_residual_test <- NA
        }
        
        if(No_Gp1>7){
            gp1_norm_res <- suppressWarnings(cvm.test(Gp1$res)$p.value)
        }else{
            gp1_norm_res <- NA
        }
        
        if(No_Gp2>7){
            gp2_norm_res <- suppressWarnings(cvm.test(Gp2$res)$p.value)
        }else{
            gp2_norm_res <- NA
        }
        
        testresults <- c(a[1], gp1_norm_res, a[2], gp2_norm_res, blups_test, 
                rotated_residual_test)
    }
    else{
        testresults <- c(a[1], NA, a[2], NA, NA, NA)
    }
    
    return(testresults)
}
##------------------------------------------------------------------------------