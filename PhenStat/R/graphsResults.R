# Q-Q Residual plots for each genotype
qqplotGenotype<-function(phenList, phenTestResult){
    if(is(phenList,"PhenList")) {
        x <- phenList$phendata     
        
    } else {
        x <- as.data.frame(phenList)
    }
    if(is(phenTestResult,"PhenTestResult")) {
        modeloutput=phenTestResult$modelOutput
    }
    else{
        modeloutput=phenTestResult
    }
    if (!is(x$Genotype)) stop("Genotype values are not defined") 
    
    res=resid(modeloutput)
    data_all= data.frame(x, res)
    a=levels(data_all$Genotype)
    par(mfrow = c(1, 2))    
    Gp1 = subset(data_all, data_all$Genotype==a[1])
    Gp2 = subset(data_all, data_all$Genotype==a[2])
    qqnorm(Gp1$res, main=" ")
    qqline(Gp1$res)
    legend("topleft", a[1], cex=1.3, bty="n")
    close.screen(n=1, all.screens = FALSE)
    qqnorm(Gp2$res, main=" ")
    qqline(Gp2$res)
    legend("topleft", a[2], cex=1.3, bty="n")
}

#-------------------------------------------------------------------------
# Predicted versus residual plots split by genotype
plotResidualPredicted<-function(phenList, phenTestResult){
    if(is(phenList,"PhenList")) {
        x <- phenList$phendata     
        
    } else {
        x <- as.data.frame(phenList)
    }
    if(is(phenTestResult,"PhenTestResult")) {
        modeloutput=phenTestResult$modelOutput
    }
    else{
       modeloutput=phenTestResult
    }
    if (!is(x$Genotype)) stop("Genotype values are not defined") 
    
    pred = predict(modeloutput)
    res=resid(modeloutput)
    data_all= data.frame(x, res, pred)
    a=levels(x$Genotype)
    genotype_no=length(a)
    par(mfrow = c(1, 2))    
    Gp1pred = subset(data_all, data_all$Genotype==a[1])
    Gp2pred = subset(data_all, data_all$Genotype==a[2])
    plot(x=Gp1pred$pred, y=Gp1pred$res, xlab="Predicted", ylab="Residuals")
    legend("topleft", a[1], cex=1.3, bty="n")
    plot(x=Gp2pred$pred, y=Gp2pred$res, xlab="Predicted", ylab="Residuals")
    legend("topleft", a[2], cex=1.3, bty="n")
}  


#-----------------------------------------------------------------------------------
# QQ plot of blups

qqplotRandomEffects<-function(phenList, phenTestResult, keep_batch=NULL){
    if(is(phenList,"PhenList")) {
        x <- phenList$phendata     
        
    } else {
        x <- as.data.frame(phenList)
    }
    if(is(phenTestResult,"PhenTestResult")) {
        modeloutput=phenTestResult$modelOutput
        if (is.null(keep_batch)) keep_batch <- phenTestResult$batchEffect
    }
    else{
        modeloutput=phenTestResult
    }
    if (is.null(keep_batch))
    stop("Please make test for the batch effect and provide TRUE/FALSE value")
    
    if(keep_batch){
        blups=ranef(modeloutput)
        qqnorm(blups[ ,1])
        qqline(blups[ ,1])

    }
    else {
        message("Diagnostics on random effects not relevant as variation between Assay.Date (batches) was not significant")
    }
}
  
#----------------------------------------------------------------------------------------------------------
# Residue versus batch split by genotype

boxplotResidualBatch<-function(phenList, phenTestResult){
    if(is(phenList,"PhenList")) {
        x <- phenList$phendata     
        
    } else {
        x <- as.data.frame(phenList)
    }
    if(is(phenTestResult,"PhenTestResult")) {
        modeloutput=phenTestResult$modelOutput
    }
    else{
        modeloutput=phenTestResult
    }
    if (!is(x$Genotype)) stop("Genotype values are not defined") 
    
    res=resid(modeloutput)
    data_all= data.frame(x, res)
    a=levels(x$Genotype)
    par(mfrow = c(1, 2))    
    Gp1 <- subset(data_all, data_all$Genotype==a[1])
    Gp2 <- subset(data_all, data_all$Genotype==a[2])
    with(Gp1, boxplot(res~Assay.Date, ylab="Residuals", xlab="Batch",names=NULL))
    legend("bottomleft", a[1], cex=1.3, bty="n")
    with(Gp2, boxplot(res~Assay.Date, ylab="Residuals", xlab="Batch",names=NULL))
    legend("bottomleft", a[2], cex=1.3, bty="n")    
}

#-------------------------------------------------------------------------------
# Q-Q plot rotated residuals

qqplotRotatedResiduals<-function(phenList, phenTestResult, keep_batch=NULL){
    if(is(phenList,"PhenList")) {
        x <- phenList$phendata     
        
    } else {
        x <- as.data.frame(phenList)
    }
    if(is(phenTestResult,"PhenTestResult")) {
        modeloutput=phenTestResult$modelOutput
        if (is.null(keep_batch)) keep_batch <- phenTestResult$batchEffect
    }
    else{
        modeloutput=phenTestResult
    }
    if (is.null(keep_batch))
    stop("Please make test for the batch effect and provide TRUE/FALSE value")
    
    x[, c("Genotype", "Gender", "Assay.Date")] = lapply(x[, c("Genotype", "Gender", "Assay.Date")], factor)
    
    if(!keep_batch){
        stop("Diagnostics on rotated residues not relevant as variation between Assay.Date (batches) was not significant")
    }else{
        sdests = exp(attr(modeloutput$apVar, "Pars"))           #extract variance estimates
        Zbat = model.matrix(~ Assay.Date, model.frame( ~ Assay.Date, modeloutput$groups))    #create random effects design matrix
        ycov = (Zbat %*% t(Zbat)) * sdests["reStruct.Assay.Date"]^2 + diag(rep(1,nrow(modeloutput$groups))) * sdests["lSigma"]^2    #create estimated cov(y)
        Lt = chol(solve(ycov))  #Cholesky decomposition of inverse of cov(y) (see Houseman '04 eq. (2))
        unrotres =  modeloutput$residuals[, "fixed"]    #unrotated residuals
        rotres = Lt %*%  modeloutput$residuals[, "fixed"]    #rotated residuals
        par(mfrow = c(1, 2))
        qqnorm(unrotres, main = "Unrotated residuals")
        qqline(unrotres)
        qqnorm(rotres, main = "Rotated residuals")
        qqline(rotres)
    }
}








