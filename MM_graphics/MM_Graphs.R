#############################################
#Raw data and Diagnostic plots 
#############################################
#test datasets
#setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBAU")
#data2=read.csv("MBAU_data.csv", header=TRUE, sep=",", dec=".")
#dataMaleonly=read.csv("MBAU_data_maleonly.csv", header=TRUE, sep=",", dec=".")


#for  diagnostic plots
#Function D1:  Raw data boxplot split by gender and genotype
boxplot_bygenotype<-function(dataset,depVariable,v_name){
	numberofgenders=length(levels(dataset$Gender))
	if(numberofgenders==2){
		Male = subset(dataset, Gender=="Male")
		Female= subset(dataset, Gender=="Female")
		par(mfrow = c(1, 2))
		boxplot(Male[ , depVariable]~Male$Genotype, ylab=v_name, xlab="Genotype")
		legend("topright", "Male", cex=1.3, bty="n")
		boxplot(Female[ , depVariable]~Female$Genotype, ylab=v_name, xlab="Genotype")
		legend("topright", "Female", cex=1.3, bty="n")
	}else{
		par(mfrow=c(1,1))
		boxplot(dataset[ ,depVariable]~dataset$Genotype, ylab=v_name, xlab="Genotype")	
	}
}	


#------------------------------------------------------------------------------
#Function D1b:    boxplot by gender, genotype and date
#Not fully tested - issues with R version and compatability between libraries.  Needs development when R version established. 
#requires R version > 2.15.3

boxplot_bygenotype_splitbybatch<-function(dataset,depVariable){
	numberofgenders=length(levels(dataset$Gender))
	require(ggplot2)
	require(grid)
		if(numberofgenders==2){
		Male = subset(dataset, Gender=="Male")
		Female= subset(dataset, Gender=="Female")
		par(mfrow = c(1, 2))
		p1=ggplot(data = Male, aes(x = Assay.Date, y = depVariable)) + geom_boxplot(aes(fill = Genotype), width = 0.8) + theme_bw()
		#legend("topright", "Male", cex=1.3, bty="n")
		p2=ggplot(data = Female, aes(x = Assay.Date, y = depVariable)) + geom_boxplot(aes(fill = Genotype), width = 0.8) + theme_bw()
		#legend("topright", "Female", cex=1.3, bty="n")
		multiplot(p1, p2, cols=2)
			}else{
		par(mfrow=c(1,1))
		ggplot(data = dataset, aes(x = Assay.Date, y = depVariable)) + geom_boxplot(aes(fill = Genotype), width = 0.8) + theme_bw()
	}
}	

#boxplot_bygenotype_splitbybatch(data2, "Fat.Mass")
#http://stackoverflow.com/questions/8320462/ggplot2-how-to-adjust-fill-colour-in-a-boxplot-and-change-legend-text
#http://stats.stackexchange.com/questions/11406/boxplot-with-respect-to-two-factors-using-ggplot2-in-r


#---------------------------------------------------------------------------
#Function D2: QQ Residual plots for each genotype for each gender
ResidualQQplot_bygenotype<-function(dataset, depVariable){
	modeloutput=finalmodel(dataset,depVariable)
	res=resid(modeloutput)
	data_all= data.frame(dataset, res)
	a=levels(data_all$Genotype)
	par(mfrow = c(1, 2))	
	Gp1 = subset(data_all, Genotype==a[1])
	Gp2 = subset(data_all, Genotype==a[2])
	qqnorm(Gp1$res, main=" ")
	qqline(Gp1$res)
	legend("topleft", a[1], cex=1.3, bty="n")
	close.screen(n=1, all.screens = FALSE)
	qqnorm(Gp2$res, main=" ")
	qqline(Gp2$res)
	legend("topleft", a[2], cex=1.3, bty="n")
}

#-------------------------------------------------------------------------
#Function D3: Predicted versus residual plots split by genotype
ResVpred_bygenotype<-function(dataset,depVariable){
	modeloutput=finalmodel(dataset,depVariable)
	pred = predict(modeloutput)
	res=resid(modeloutput)
	data_All= data.frame(dataset, res, pred)
	a=levels(dataset$Genotype)
	genotype_no=length(a)
	par(mfrow = c(1, 2))	
	Gp1pred = subset(data_All, Genotype==a[1])
	Gp2pred = subset(data_All, Genotype==a[2])
	plot(x=Gp1pred$pred, y=Gp1pred$res, xlab="Predicted", ylab="Residual")
	legend("topleft", a[1], cex=1.3, bty="n")
	plot(x=Gp2pred$pred, y=Gp2pred$res, xlab="Predicted", ylab="Residual")
	legend("topleft", a[2], cex=1.3, bty="n")
}	

#-----------------------------------------------------------------------------------
#Function D4:  Body weight versus dependent variable scatter plot

weight_versus_depVariable<-function(dataset, depVariable){
	require(car)
	model.formula <- as.formula(paste(depVariable, "~", paste("Weight", "Genotype", sep= "|")))
	scatterplot(data=dataset, model.formula)
}

#----------------------------------------------------------------------------------------------------
#Function D5:    QQ plot of Blups
		
blups_plot<-function(dataset, depVariable){
	modeloutput=finalmodel(dataset,depVariable)
	blups=ranef(modeloutput)
	qqnorm(blups[ ,1])
	qqline(blups[ ,1])
}

#----------------------------------------------------------------------------------------------------------
#Function D6: Residue versus batch split by genotype
		
ResidualVersusBatch_bygenotype<-function(dataset, depVariable){
	modeloutput=finalmodel(dataset,depVariable)
	res=resid(modeloutput)
	data_all= data.frame(dataset, res)
	a=levels(dataset$Genotype)
	par(mfrow = c(1, 2))	
	Gp1 <- subset(data_all, Genotype==a[1])
	Gp2 <- subset(data_all, Genotype==a[2])
	with(Gp1, boxplot(res~Assay.Date, ylab="residues", xlab="Batch"))
	legend("bottomleft", a[1], cex=1.3, bty="n")
	with(Gp2, boxplot(res~Assay.Date, ylab="residues", xlab="Batch"))
	legend("bottomleft", a[2], cex=1.3, bty="n")	
}

#-------------------------------------------------------------------------------
#Function D7: Rotated residuals QQ plot

rot_residuals<-function(dataset, depVariable){
	dataset[, c("Genotype", "Gender", "Assay.Date")] = lapply(dataset[, c("Genotype", "Gender", "Assay.Date")], factor)
	keep_batch=testing_batch(dataset, depVariable)
	if(!keep_batch){
		stop("Diagnostics on rotated residues not relevant as variation between Assay.Date was not significant")
	}else{
		modeloutput=finalmodel(dataset,depVariable)             #fit lmm
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


