
#-------------------------------------------------------------------------------------------------------------------------
#Function EnsureWTisRefLevel  is a function designed to ensure when one of the genotypes is wildtype then wildtype will always be the reference level and this will ensure the comparison is always knockout relative to wildtype.  
#In the situation where there is no wildtype data then it doesn't matter which is the reference.  In the MM output group1 genotype will always be the reference set and thus it can be determined from output.
#Placed in the point where the output is captured as it doesn't effect any of the steps up until that point.  

EnsureWTisRefLevel<-function(dataset){
	Genotype_levels=levels(dataset$Genotype)
	
	if(sum(grepl("+/+", Genotype_levels, fixed=TRUE))==1){
		dataset$Genotype=relevel(dataset$Genotype, ref="+/+")
		return(dataset)
	}
}


#-----------------------------------------------------
#Function testing_batch: function to test batch
# A full loaded model formula (ie one with all variables of interest) is assembled.  Weight cannot be in the model as a dependent variable if it is the independent variable - hence the if else questions to build the model.
# Then the model formula are used to build two models one with batch included as a random effect and hence uses a mixed model and one where batch is not included and uses linear regression
# The two models are compared via anova to test the null hypothesis that batch is not significant (Ho) where the alternative hypothesis (Ha) is that batch is significant.  If the p value <0.05,  then the we reject the null hypothesis and accept the alternative and say we need to keep batch
# the p-value for batch test is divided by 2 as the test is on the boundary of the parameter space for a variance and the null distribution of the likelihood ratio test statistics follows a mixture of chi-squared distributions with equal weight of 0.5.
# As the output is yes or no on keep_batch then if p<0.05,  we say batch is significant  and a mixed model approach should be used in subsequent work.

testing_batch<-function(dataset, depVariable){
	require(nlme)
	model.formula <- model_Formula(dataset, depVariable)
	model_MM = do.call("lme", args = list(model.formula, random=~1|Assay.Date, data = dataset, na.action="na.omit", method="REML"))
	model_withoutbatch <- do.call("gls", args=list(model.formula, dataset, na.action="na.omit"))
	p.value.batch <-(anova(model_MM, model_withoutbatch)$p[2])/2
	keep_batch= p.value.batch<0.05
	return(keep_batch)
}

#---------------------------------------------------
#Function testing_variance: function to test equal variance

# A full loaded model forumula(ie one with all variables of interest) is assembled.  Weight cannot be in the model as a dependent variable if it is the independent variable - hence the if else questions to build the model.
# Then the models formula are used to build two models one with homogenous variance and one with heterogenous variance
# The two models are compared via anova to test the null hypothesis that variance is equal (Ho) where the alternative hypothesis (Ha) is that variance is not equal.  If the p value <0,05,  the null hypotehsis is rejected and the alternative accepted.
# As the output is yes or no on keep_equalvar then if p>0.05,  we say variance should be equal

testing_variance<-function(dataset, depVariable){
	require(nlme)
	model.formula <- model_Formula(dataset, depVariable)
	model_equalvariance=do.call("lme", args=list(model.formula, random=~1|Assay.Date, dataset, na.action="na.omit", method="REML") )
	model_hetvariance= do.call("lme", args=list(model.formula, random=~1|Assay.Date, dataset, weights=varIdent(form=~1|Genotype), na.action="na.omit", method="REML"))
	p.value=(anova(model_equalvariance, model_hetvariance)$p[2])
	keep_equalvar= p.value>0.05
	return(keep_equalvar)
}

#------------------------------------------------------------------
# Function model_forFIXEDtest: function to build model following variance and batch test for testing of fixed effects
# The keep_batch and keep_equal variance functions are called to build the correct models via a series of if else rules.

model_forFIXEDtest<-function(dataset, depVariable){
	require(nlme)
	#test whether batch is significant
	keep_batch=testing_batch(dataset,depVariable)
	#test whether variance is significant
	keep_equalvar=testing_variance(dataset, depVariable)
	# outputs the fully loaded model for analysis
	model.formula <- model_Formula(dataset, depVariable)
	#build correct model depending on the output of the two questions (batch and variance)
	
	if(keep_batch && keep_equalvar){
		return(model = lme(model.formula, random=~1|Assay.Date, dataset, na.action="na.omit", method="REML"))
	}else if(keep_batch && !keep_equalvar){
		return(model=lme(model.formula, random=~1|Assay.Date, dataset, weights=varIdent(form=~1|Genotype), na.action="na.omit", method="REML"))
	}else if(!keep_batch && keep_equalvar){
		return(model=gls(model.formula, dataset, na.action="na.omit"))
	}else if(!keep_batch && !keep_equalvar){
		return(model=gls(model.formula, weights=varIdent(form=~1|Genotype), dataset, na.action="na.omit"))
	}
}

#----------------------------------------------------------------------------------
#Function testing_genotype_effect: The goal of the function is to give a pvalue of whether genotype effect is significant by compare the genotype model with the null model with the anova function
# The model_formula_null and model_formula_genotype are called to define the models for comparison.
# The keep_batch and keep_equal variance function are called this allows if but rules to build the correct model comparison ie if batch is included then we use a mixed model etc
# For each possible combination,  then an anova model is used to report the pvalue.
# For testing the genotype effect we use method= ML


testing_genotype_effect<-function(dataset, depVariable){
	require(nlme)
	#following needed to determine model general structure
	keep_batch=testing_batch(dataset,depVariable)
	keep_equalvar=testing_variance(dataset, depVariable)
	
	#call functions to determine model.formula
	model_formula_null = final_null_model(dataset,depVariable)
	model_formula_genotype = final_genotype_model(dataset, depVariable)
	
	
	if(keep_batch && keep_equalvar){
		model_genotype=do.call("lme", args = list(model_formula_genotype, random=~1|Assay.Date, dataset, na.action="na.omit", method="ML"))
		model_null=do.call("lme", args=list(model_formula_null, dataset,random=~1|Assay.Date, na.action="na.omit",  method="ML"))
		return(p.value=(anova(model_genotype, model_null)$p[2]))
	}else if(keep_batch && !keep_equalvar){
		model_genotype=do.call("lme", args = list(model_formula_genotype, random=~1|Assay.Date, dataset,weights=varIdent(form=~1|Genotype), na.action="na.omit", method="ML"))
		model_null=do.call("lme", args=list(model_formula_null, dataset, random=~1|Assay.Date,weights=varIdent(form=~1|Genotype), na.action="na.omit",  method="ML"))
		return(p.value=(anova(model_genotype, model_null)$p[2]))
	}else if(!keep_batch && !keep_equalvar){
		model_genotype=do.call("gls", args = list(model_formula_genotype,  dataset,weights=varIdent(form=~1|Genotype),method="ML", na.action="na.omit"))
		model_null=do.call("gls", args=list(model_formula_null, dataset,weights=varIdent(form=~1|Genotype),method="ML", na.action="na.omit"))
		return(p.value=(anova(model_genotype, model_null)$p[2]))
	}else if(!keep_batch && keep_equalvar){
		model_genotype=do.call("gls", args = list(model_formula_genotype,  dataset, method="ML", na.action="na.omit"))
		model_null=do.call("gls", args=list(model_formula_null, dataset,method="ML", na.action="na.omit"))
		return(p.value=(anova(model_genotype, model_null)$p[2]))
	}
}

#-----------------------------------------------------------------------------------------
#Function finalmodel: following function returns the final model which is needed for the diagnostic plots
#Goal of this function is to return a model output that is the final model following all the previous tests which can be used to generate diagnostics plots and output the final model details.
# The model_formula_null and model_formula_genotype are called to define the models for comparison.
# The keep_batch and keep_equal variance function are called this allows if but rules to build the correct model comparison ie if batch is included then we use a mixed model etc
# For each possible combination,  then the model is fitted to the data and the model reported as the output.
#this is a separate function to above as for the estimates of the fixed effects we use method=REML whilst for the genotype test we used method= ML
#na.action = na.exclude as in this form the residue calculated from the model output have the same length as the original datafile

finalmodel<-function(dataset, depVariable){
	
	require(nlme)
	#following needed to determine model general structure
	keep_batch=testing_batch(dataset,depVariable)
	keep_equalvar=testing_variance(dataset, depVariable)
	
	#call functions to determine model.formula
	model_formula_null = final_null_model(dataset,depVariable)
	model_formula_genotype = final_genotype_model(dataset, depVariable)
	
	if(keep_batch && keep_equalvar){
		return(model_genotype=do.call("lme", args = list(model_formula_genotype, random=~1|Assay.Date, dataset, na.action="na.exclude", method="REML")))
	}else if(keep_batch && !keep_equalvar){
		return(model_genotype=do.call("lme", args = list(model_formula_genotype, random=~1|Assay.Date, dataset,weights=varIdent(form=~1|Genotype), na.action="na.exclude", method="REML")))
	}else if(!keep_batch && !keep_equalvar){
		return(model_genotype=do.call("gls", args = list(model_formula_genotype,  dataset,weights=varIdent(form=~1|Genotype), na.action="na.exclude")))
	}else if(!keep_batch && keep_equalvar){
		return(model_genotype=do.call("gls", args = list(model_formula_genotype,  dataset, na.action="na.exclude")))
	}
}

#---------------------------------------



countDataPoints<-function(dataset, depVariable){
	print(sum(is.finite(dataset[ , depVariable])))
}


