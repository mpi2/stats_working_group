# TODO: Add comment
# 
# Author: nk3
###############################################################################
if  (!library(partitions, logical.return = TRUE)) {	install.packages("partitions") }

if  (!library(Epi, logical.return = TRUE)) { install.packages("Epi") }

if  (!library(epitools, logical.return = TRUE)) {	install.packages("epitools") }

require(logistf)
require(epitools)
require(partitions)
require(Epi)

#Function to return formula for model fitting for the logistif regression
#Two arguments 
#First: modelType can be "full" or "null" 
#Second: depVariableString is a string to represent the dependent Variable eg "Sodium"
#Example use model_Formula("full", "Sodium")

model_Formula <- function(modelType, depVariableString){
	
	if(modelType=="full"){
		model.formula <- as.formula(paste(depVariableString,  "Sex", sep="~"))
	}else if(modelType=="null"){ 
		model.formula <- as.formula(paste(depVariableString, "~", "1",  sep= " ")) 
	}
	return(model.formula)
}


#Function to  return the stage two testing (pvalue, estimated coefficient and SE on estimated coefficient).  
#Here we test the differences in abnormality rate between the male and female knockout data for a dependent variable using logistf
#Assumes the variable has been recoded appropriately to 0 and 1
#Assumes the dataframe has been prepared and only knockout data is presented
#Assumes the dataframe has a column called Sex and a column for the dependentVarible labelled with a string depVariableString

LR_KOdata<-function(KOdata, depVariableString){
	
	Signal2=length(levels(as.factor(KOdata[ ,depVariableString]))) ##assessing whether signal within the knockout mice to allow model fitting
	
	if(Signal2==2){
		
		formula_full_stage2=model_Formula(modelType="full", depVariableString)  #returns the test model
		model_full_stage2 <- do.call("logistf",	args=list(formula_full_stage2, data=KOdata, na.action="na.omit"))  #fits the test model to the knockoutdata
		LR_KO_pval=logistftest(model_full_stage2)$prob   # strategy needed to test effect  as you cannot specify an intercept only model with logistf you cannot anova(test, null)
		finalmodel=summary(model_full_stage2)				
		output=c(LR_KO_pval, finalmodel$coefficients[2], format(sqrt(diag(vcov(finalmodel))),scientific=FALSE)[2])
		names(output)=c("LR_KO_pval", "LR_KO_est_coeff", "LR_KO_coeff_SE")
	}else{
		output=c(1, 0,0)  #when there is no abnormality signal returns a p value of 1 and 0 for estimate and SE
		names(output)=c("LR_KO_pval", "LR_KO_est_coeff", "LR_KO_coeff_SE")
	}		
	return(output)	
}



###LR_KOalpha star  p value which is the most extreme value possible arising as a function of the number of abnormal calls and number of readings taken. 

##Using the annotation used within this code:
#' Females:             Males: 
#'   
#' |  |KO  |WT  |    |  |  |KO  |WT  |    |
#' |--|:--:|:--:|:--:|  |--|:--:|:--:|:--:|
#' |1 |y_f |    |yx_f|  |1 |y_m |    |yx_m|
#' |0 |    |    |    |  |0 |    |    |    |
#' |  |n2_f|n1_f|    |  |  |n2_m|n1_m|    |

##The 2 by 2 table for stage 2 would be: 
#' |  |KO _M |KO_F  |  
#' |--|:-  -:|:--:  |
#' |1 |y_m   |  y_f |  
#' |0 |      |      |      
#' |  |n2_m  |n2_f  |    


LR_KOalphaStar<-function(y_m, y_f, n2_m, n2_f){
	
	#build a dataframe based on the most extreme values being loaded into the males 	
	if((y_m+y_f)<=n2_m){  			#number of abnormalities is less than the number of male mice
		
		#' |  |KO _M                |KO_F     |  
		#' |--|                     |         |
		#' |1 |(y_m +y_f )          |         | 
		#' |0 |(n2_m)-(y_m +y_f )   |n2_f     |    		
		
		male_df <- data.frame(AbnormalityCall=c(rep(x=0, times=((n2_m)-(y_m +y_f ))), rep(x=1, times=(y_m +y_f ))), Sex=rep(x="Male", times=n2_m))
		female_df <- data.frame(AbnormalityCall=rep(x=0, times=n2_f), Sex=rep(x="Female", times=n2_f))
		male_extreme=rbind(male_df, female_df)
		
	}else { #number of abnormalities is greater than the number of male mice
		#if	((y_m+y_f)>n2_m)	
		#' |  |KO _M          |KO_F                    |  
		#' |--|               |                        |
		#' |1 |n2_m	          |(y_m+y_f)-n2_m          |
		#' |0 |		   	   	|n2_f- ((y_m+y_f)-n2_m )   |  
		
		male_df <- data.frame(AbnormalityCall=rep(x=1, times=n2_m), Sex=rep(x="Male", times=n2_m))
		female_df <- data.frame(AbnormalityCall=c(rep(x=1, times=((y_m+y_f)-n2_m)), rep(x=0, times=(n2_f- ((y_m+y_f)-n2_m )))), Sex=rep(x="Female", times=n2_f))
		male_extreme=rbind(male_df, female_df)	
	}		
	require(logistf)
	
	Signal2=length(levels(as.factor(male_extreme[ ,"AbnormalityCall"]))) ##assessing whether sufficient signal within the male_extreme table to allow model fitting
	
	if(Signal2==2){
		
		formula_full_stage2=model_Formula(modelType="full", depVariableString="AbnormalityCall")  #returns the model formula Y~Sex
		model_full_stage2 <- do.call("logistf",	args=list(formula_full_stage2, data=male_extreme, na.action="na.omit"))  #fits the test model to the male extreme data
		Male_extreme_pval=logistftest(model_full_stage2)$prob   # strategy needed to test effect  as you cannot specify an intercept only model with logistf you cannot anova(test, null)
		
	}else{
		Male_extreme_pval=1  #when there is no abnormality signal returns a p value of 1 	
	}		
	#build a dataframe based on the most extreme values being loaded into the females 	
	if((y_m+y_f)<=n2_f){  			#number of abnormalities is less than the number of female mice
		
		#' |  |KO _M         |KO_F              |  
		#' |--|              |                  |
		#' |1 |              |(y_m +y_f )       | 
		#' |0 |(n2_m)        |n2_f-(y_m +y_f )  |     
		female_df <- data.frame(AbnormalityCall =c(rep(x=0, times=(n2_f-(y_m +y_f ))), rep(x=1, times=(y_m +y_f ))), Sex=rep(x="Female", times=n2_f))
		male_df <- data.frame(AbnormalityCall =rep(x=0, times=n2_m), Sex=rep(x="Male", times=n2_m))
		female_extreme=rbind(male_df, female_df)
		
		
	}else { #number of abnormalities is greater than the number of female mice
		#if	((y_m+y_f)>n2_f)	
		
		#' |  |KO _M        		  |KO_F                    |  
		#' |--|              		  |                        |
		#' |1 |(y_m+y_f)-n2_f 	      |n2_f                    |
		#' |0 |	n2_m-((y_m+y_f)-n2_f) |	                       |  
		
		female_df <- data.frame(AbnormalityCall=rep(x=1, times=n2_f), Sex=rep(x="Female", times=n2_m))
		male_df <- data.frame(AbnormalityCall=c(rep(x=1, times=((y_m+y_f)-n2_m)), rep(x=0, times=(n2_m-((y_m+y_f)-n2_f )))), Sex=rep(x="Male", times=n2_f))
		female_extreme=rbind(male_df, female_df)	
	}
	
	Signal2=length(levels(as.factor(female_extreme[ ,"AbnormalityCall"]))) ##assessing whether sufficient signal within the female_extreme table to allow model fitting
	if(Signal2==2){	
		formula_full_stage2=model_Formula(modelType="full", depVariableString="AbnormalityCall")  #returns the model formula Y~Sex
		model_full_stage2 <- do.call("logistf",	args=list(formula_full_stage2, data=female_extreme, na.action="na.omit"))  #fits the test model to the female extreme data
		Female_extreme_pval=logistftest(model_full_stage2)$prob   # strategy needed to test effect  as you cannot specify an intercept only model with logistf you cannot anova(test, null)	
	}else{
		Female_extreme_pval=1  #when there is no abnormality signal returns a p value of 1 	
	}		
	output=min(Male_extreme_pval, Female_extreme_pval)
	#names(output) = "LR_KO_alphaStar"
	return(output)	
}






#' Notations
#'  
#' Females:             Males: 
#'   
#' |  |KO  |WT  |    |  |  |KO  |WT  |    |
#' |--|:--:|:--:|:--:|  |--|:--:|:--:|:--:|
#' |1 |y_f |    |yx_f|  |1 |y_m |    |yx_m|
#' |0 |    |    |    |  |0 |    |    |    |
#' |  |n2_f|n1_f|    |  |  |n2_m|n1_m|    |



#' Cochran-Mantel-Haenszel Test
#'
#' Tests the proprtion of rare events difference between the KO and WT groups, stratified by Gender.
#' The test is one-sided. Tests whether the KO group has larger proportion than the WT.
#' 
#' @param y_f, yx_f, n2_f, n1_f, y_m, yx_m, n2_m, n1_m
#' Margins and values of the 2x2x2 contingancy table (see notations above)
#' 
#' @param example
#' Logical. If TRUE the example from [*] is analysed (useful for debugging or demonstration).
#'
#' @return Numeric vector with 3 values :alpha-star, p-value and mid p-value.
#' 
#' @example 
#' MH_test(example=T)
#' 
#' [*] book: Nonparametric Statistical Methods, 3e Hollander, Wolfe & Chicken, chapter 10) 

MH_test <- function (y_f, yx_f, n2_f, n1_f, y_m, yx_m, n2_m, n1_m, example = F, alpha_star = "old")
{
	if (example)
	{
		z <- array(c(2, 1, 2, 5, 1, 5, 4, 1), dim=c(2, 2, 2))
		y_f  <- z[1,1,1]
		n2_f <- sum(z[,1,1])
		n1_f <- sum(z[,2,1])
		yx_f <- sum(z[1,,1])
		y_m  <- z[1,1,2]
		n2_m <- sum(z[,1,2])
		n1_m <- sum(z[,2,2])
		yx_m <- sum(z[1,,2])
	}
	else
		z <- array(c(y_f,n2_f-y_f,yx_f-y_f,n1_f-(yx_f-y_f),y_m,n2_m-y_m,yx_m-y_m,n1_m-(yx_m-y_m)), dim=c(2, 2, 2))
	
	# alpha star:
	# high limit of the support
	hi_vec <- c(min(n2_f,yx_f),min(n2_m,yx_m))
	
	if (alpha_star == "new")
	{
		z_as <- array(c(n2_f ,0 ,max(0,yx_f-n2_f) ,n1_f-max(0,yx_f-n2_f)  ,
						n2_m ,0 ,max(0,yx_m-n2_m) ,n1_m-max(0,yx_m-n2_m)  ), dim=c(2, 2, 2))
	}
	else
	{
		z_as  <- array(c(hi_vec[1],n2_f-hi_vec[1],yx_f-hi_vec[1],n1_f-(yx_f-hi_vec[1]),
						hi_vec[2],n2_m-hi_vec[2],yx_m-hi_vec[2],n1_m-(yx_m-hi_vec[2])), dim=c(2, 2, 2))
	}
	as <- mantelhaen.test(x = z_as,alternative = "greater", exact = T)$p
	
	# p-value:
	z <- array(c(y_f,n2_f-y_f,yx_f-y_f,n1_f-(yx_f-y_f),y_m,n2_m-y_m,yx_m-y_m,n1_m-(yx_m-y_m)), dim=c(2, 2, 2))
	out <- mantelhaen.test(x = z,alternative = "greater", exact = T)
	pv <- out$p
	
	# mid-p-value:
	s <- out$s
	if (s==sum(hi_vec))
		mid_pv <- pv/2
	else
	{
		z_1 <- z
		if (z[1,1,1]<hi_vec[1])
		{
			z_1[1,1,1] <- z[1,1,1] + 1
			z_1[2,1,1] <- z[2,1,1] - 1
			z_1[1,2,1] <- z[1,2,1] - 1
			z_1[2,2,1] <- z[2,2,1] + 1
		}
		else
		{
			z_1[1,1,2] <- z[1,1,2] + 1
			z_1[2,1,2] <- z[2,1,2] - 1
			z_1[1,2,2] <- z[1,2,2] - 1
			z_1[2,2,2] <- z[2,2,2] + 1
		}
		pv_1 <- mantelhaen.test(x = z_1,alternative = "greater", exact = T)$p
		mid_pv <- (pv+pv_1)/2
	}
	out <- c(as,pv,mid_pv)
	names(out) <- c("alpha_star","pv","mid_pv")
	return(out)
}



	
	
#' Zelen test
#'
#' Tests whether the difference in proprtion (KO vs. WT) differ across gender (interaction).
#' The test is two sided.
#' 
#' @param y_f, yx_f, n2_f, n1_f, y_m, yx_m, n2_m, n1_m
#' Margins and values of the 2x2x2 contingancy table (see notations above)
#' 
#' @param example
#' Logical. If TRUE the example from [*] is analysed (useful for debugging or demonstration).
#'
#' @return Numeric vector with 3 values :alpha-star, p-value and mid p-value.
#' 
#' @example 
#' zelen_test(example = T)
#' 
#' Base on zelen.test() function by Eric Chicken from package NSM3.
	
	
	zelen_test <- function (y_f, yx_f, n2_f, n1_f, y_m, yx_m, n2_m, n1_m, example=F) 
	{
		# Zelen's test.
		# Based on chapter 10 of:
		#
		#   Nonparametric Statistical Methods, 3e
		#   Hollander, Wolfe & Chicken 
		#
		# The data z is an array of k 2x2 matrices.
		# Small data sets only!
		# Uses package "partitions".
		#
		# Inefficiently programmed by Eric Chicken, November 2012.
		
		if(example)
			z <- array(c(2, 1, 2, 5, 1, 5, 4, 1), dim=c(2, 2, 2))
		else
			z <- array(c(y_f,n2_f-y_f,yx_f-y_f,n1_f-(yx_f-y_f),y_m,n2_m-y_m,yx_m-y_m,n1_m-(yx_m-y_m)), dim=c(2, 2, 2))
		
		s <- sum(z[1, 1, ])
		k <- dim(z)[3]
		
		# blockparts is from package "partitions".  This is where large data
		# sets will be an issue.
		# Make sure that each part of the sum is no more than the column or
		# row margin total.
		bp <- numeric(0)
		for(i in 1:k) bp <- c(bp, min(sum(z[1,,i]),sum(z[,1,i])))
		a <- blockparts(bp, s)
		
		y <- numeric(0)
		for(i in 1:dim(a)[2])
		{
			is.tau.0 <- T
			x <- numeric(0)
			for(j in 1:k)
			{
				O.11 <- a[j, i]
				O.12 <- sum(z[1, , j]) - O.11
				O.21 <- sum(z[, 1, j]) - O.11
				O.22 <- sum(z[2, , j]) - O.21
				tau <- matrix(c(O.11, O.12, O.21, O.22), nrow=2, byrow=T)
				if(sum(tau == z[, , j]) < 4) is.tau.0 <- F
				n1 <- O.11 + O.12
				n2 <- O.21 + O.22
				n.1 <- O.11 + O.21
				n <- n1 + n2
				x.j <- choose(n1, O.11) * choose(n2, O.21) / choose(n, n.1)
				x <- c(x, x.j)
			}
			if(is.tau.0) tau.0 <- i
			y <- c(y, prod(x))
		}
		y <- y / sum(y)
		pv <- sum(y[y<=y[tau.0]])
		#   cat("\n")
		#   cat("Zelen's test:") 
		#   cat("\n")
		#   cat(paste("P = ", round(p, r), sep=""))
		#   cat("\n")
		alpha_star <- min(y)
		f_n <- y[y==y[tau.0]][1]
		mid_pv <- pv - 0.5*f_n
		
		out <- c(alpha_star,pv,mid_pv)
		names(out) <- c("alpha_star","pv","mid_pv")
		
		return(out)
	}
	
	


#Inputs
#df  a dataframe of knockout and wildtype data with a column for Sex, includes multiple columns for the variables of interest which have been recoded to 0/1, a Genotype column (which specifies wildtype by +/+),
#a colonyPrefix column which has a string which specifies which knockoutdata should be grouped for a gene and a column Zygosity which has a HOM/HET/WT string which specifies the version of the gene.
#requires df to have a column Assay.Date
#requires df to have two sexes
#returns a vector of 25 elements

#Version with tryCatch management of errors	
Cat_GenotypeSDeffect_pipeline2<-function(df,dependentVariable, ZygosityToTest, refGenotype="WT"){
		
			test=PhenList(df, testGenotype=ZygosityToTest, refGenotype=refGenotype, dataset.clean=TRUE, dataset.colname.batch="Assay.Date",  dataset.colname.genotype="Zygosity", outputMessages=FALSE)
			
				
				tryCatch(
						{				
							# runs the PhenStat FE method
							result=testDataset(phenList=test, depVariable=dependentVariable, outputMessages=FALSE, pThreshold=0.05, method="FE", transformValues=FALSE)
							#FEresults= vectorOutput(result)
							#as PhenStat has a standard vectorOutput format for all methods there are things we don't need when only focused on #FE analysis
							#FEresultsKeep=c(FEresults[c(1,2, 16,18,26,28,29,31,32)], as.character(Ds_Colony), as.character(Ds_Genotype), screen)
							#names(FEresultsKeep)=c("Method", "DepVariable", "Gp1Genotype", "Gp2Genotype", "FE_WTKO_Fem_ES", "FE_WTKO_Fem_pvalue", "FE_WTKO_male_ES", "FE_WTKO_male_pvalue", "PhenStat_FE_classification")
							
							#pull out the count matrix from PhenStat results
							allCountMatrices <- getCountMatrices(result)
							#They are organised as matrix as below:
							
							#' Females:             Males: 
							#'   
							#' |  |WT    |KO    |    |  |WT  |KO  |   
							#' |--|:--   |:--:  |:--:|  |:--:|:--:|
							#' |0 |[1,1] |[1,2] |    |0 |    |    |
							#' |1 |[2,1] |[2,2] |    |1 |    |    |   
							
							
							#Newcombe 95% CI for difference of two independent proportions 
							#http://www.inside-r.org/packages/cran/epi/docs/ci.pd  
							#Here the ES and CI is claculated for each table (male KO-WT and female KO-WT) separately.  To give a summary ES for stage 1, which is a calculating the main effect across both tables.  The average ES is calcualated and a CI selected by taking the min of the lower CI and the max of the upper CI from the two separate calculations.  This would be a conservative estimate. 
							
							CIoutput_Male=ci.pd(aa=allCountMatrices$male[2,1], bb=allCountMatrices$male[2,2], cc=allCountMatrices$male[1,1] ,dd=allCountMatrices$male[1,2],		method = "Nc",		alpha = 0.05, conf.level=0.95,		digits = 3,		print = FALSE,		detail.labs = FALSE )
							CIoutput_Female=ci.pd(aa=allCountMatrices$female[2,1], bb=allCountMatrices$female[2,2], cc=allCountMatrices$female[1,1] ,dd=allCountMatrices$female[1,2],		method = "Nc",		alpha = 0.05, conf.level=0.95,		digits = 3,		print = FALSE,		detail.labs = FALSE )
							CIoutput=c(CIoutput_Male[c(2, 4,5,6,7)],CIoutput_Female[c(2, 4,5,6,7)], (CIoutput_Male[5]+CIoutput_Female[5])/2, min(CIoutput_Male[6],CIoutput_Female[6]), max(CIoutput_Male[7],CIoutput_Female[7]) )   
							names(CIoutput)=c("Male_AbRateWT", "Male_AbRateKO", "Male_DiffAbRate", "Male_lower95CI", "Male_upper95CI","Fem_AbRateWT", "Fem_AbRateKO", "Fem_DiffAbRate", "Fem_lower95CI", "Fem_upper95CI" , "Stage1_AvDiffAbnRate", "Stage1_lowerCI", "Stage1_upperCI")
							
							#print(CIoutput)
							
							#stage 2 testing: Logistf regression comparing abnormality rates acrosses the sexes in knockout data only
							KOdatasetOnly=subset(result@analysedDataset, result@analysedDataset$Genotype!="WT")  #prepare knockout data for stage 2 function
							Stage2_InteractionTest=LR_KOdata(KOdatasetOnly, depVariableString=dependentVariable)
							
							#stage 2 testing: LR_KO alpha star p value
							LR_KOalphaStar=LR_KOalphaStar(y_m=allCountMatrices$male[2,2], y_f=allCountMatrices$female[2,2], n2_m=(allCountMatrices$male[2,2]+allCountMatrices$male[1,2]), n2_f=(allCountMatrices$female[2,2]+allCountMatrices$female[1,2]))
							names(LR_KOalphaStar)	=c("LR_KOalphaStar")					
							
							
							#stage 1 testing: MH_test
							Stage1_MH=MH_test(y_f=allCountMatrices$female[2,2], yx_f=(allCountMatrices$female[2,2]+allCountMatrices$female[2,1]), n2_f=(allCountMatrices$female[2,2]+allCountMatrices$female[1,2]), n1_f=(allCountMatrices$female[1,1]+allCountMatrices$female[2,1]), y_m =allCountMatrices$male[2,2], yx_m=(allCountMatrices$male[2,2]+allCountMatrices$male[2,1]), n2_m=(allCountMatrices$male[2,2]+allCountMatrices$male[1,2]), n1_m =(allCountMatrices$male[1,1]+allCountMatrices$male[2,1]), example = F)			
							names(Stage1_MH)=c("MH_alpha_star", "MH_pv", "Stage1_MH_midpv")			
							
							
							#stage 2 testing: zelen_test with one added rare event for females WT and one for males WT
							#' Females:             Males: 
							#'   
							#' |  |KO  |WT  |    |  |  |KO  |WT  |    |
							#' |--|:--:|:--:|:--:|  |--|:--:|:--:|:--:|
							#' |1 |y_f |    |yx_f|  |1 |y_m |    |yx_m|
							#' |0 |    |    |    |  |0 |    |    |    |
							#' |  |n2_f|n1_f|    |  |  |n2_m|n1_m|    |

							
							Stage2_zelen=zelen_test(y_f=allCountMatrices$female[2,2], yx_f=(allCountMatrices$female[2,2]+allCountMatrices$female[2,1]) + 1, n2_f=(allCountMatrices$female[2,2]+allCountMatrices$female[1,2]), n1_f=(allCountMatrices$female[1,1]+allCountMatrices$female[2,1]), y_m =allCountMatrices$male[2,2], yx_m=(allCountMatrices$male[2,2]+allCountMatrices$male[2,1]) + 1, n2_m=(allCountMatrices$male[2,2]+allCountMatrices$male[1,2]), n1_m =(allCountMatrices$male[1,1]+allCountMatrices$male[2,1]), example = F)			
							names(Stage2_zelen)=c("zelen_alpha_star", "zelen_pv", "zelen_midpv")	
							
							
							#Alternate stage2 testing: FE comparing abnormality rate acrosses the sexes in knockout data only
							#challenge.df  <-matrix(c(allCountMatrices$male[2,2],allCountMatrices$male[1,2],allCountMatrices$female[2,2],allCountMatrices$female[1,2]),
							#		nrow = 2,	dimnames = list(phenotype = c("1", "0"),	Sex = c("Male", "Female")))
							#Stage2_FE=fisher.test(challenge.df)$p.value
							#names(Stage2_FE)=c("FE_Stage2_KOonly")
							
							
							#CI for stage 2 test
							CI_stage2=ci.pd(aa=allCountMatrices$female[2,2], bb=allCountMatrices$male[2,2], cc=allCountMatrices$female[1,2] ,dd=allCountMatrices$male[1,2],	method = "Nc",alpha = 0.05, conf.level=0.95,		digits = 3,		print = FALSE,		detail.labs = FALSE )
							CI_stage2_output=c(CI_stage2[c(5,6,7)])   
							names(CI_stage2_output)=c("Stage2DiffAbRate", "Stage2lower95CI", "Stage2upper95CI" )

							
							CombinedResults=c(Stage1_MH, CIoutput,Stage2_InteractionTest,  CI_stage2_output, Stage2_zelen, LR_KOalphaStar)
							#print(CombinedResults)
							
						},					
						error=function(e){
							
							filename=paste(paste(dependentVariable,ZygosityToTest, "Failure",  sep="_"), "csv", sep=".")
							write.csv(x=df, file=filename)
						})	
			return(as.numeric(CombinedResults))
}
	
#Version without tryCatch management of errors	for jeremy mason
Cat_GenotypeSDeffect_pipeline2v2<-function(df,dependentVariable, ZygosityToTest, refGenotype="WT"){
	
	test=PhenList(df, testGenotype=ZygosityToTest, refGenotype=refGenotype, dataset.clean=TRUE, dataset.colname.batch="Assay.Date",  dataset.colname.genotype="Genotype", outputMessages=FALSE)
		
				# runs the PhenStat FE method
				result=testDataset(phenList=test, depVariable=dependentVariable, outputMessages=FALSE, pThreshold=0.05, method="FE", transformValues=FALSE)
				#FEresults= vectorOutput(result)
				#as PhenStat has a standard vectorOutput format for all methods there are things we don't need when only focused on #FE analysis
				#FEresultsKeep=c(FEresults[c(1,2, 16,18,26,28,29,31,32)], as.character(Ds_Colony), as.character(Ds_Genotype), screen)
				#names(FEresultsKeep)=c("Method", "DepVariable", "Gp1Genotype", "Gp2Genotype", "FE_WTKO_Fem_ES", "FE_WTKO_Fem_pvalue", "FE_WTKO_male_ES", "FE_WTKO_male_pvalue", "PhenStat_FE_classification")
				
				#pull out the count matrix from PhenStat results
				allCountMatrices <- getCountMatrices(result)
				#They are organised as matrix as below:
				
				#' Females:             Males: 
				#'   
				#' |  |WT    |KO    |    |  |WT  |KO  |   
				#' |--|:--   |:--:  |:--:|  |:--:|:--:|
				#' |0 |[1,1] |[1,2] |    |0 |    |    |
				#' |1 |[2,1] |[2,2] |    |1 |    |    |   
				
				
				#Newcombe 95% CI for difference of two independent proportions 
				#http://www.inside-r.org/packages/cran/epi/docs/ci.pd  
				#Here the ES and CI is claculated for each table (male KO-WT and female KO-WT) separately.  To give a summary ES for stage 1, which is a calculating the main effect across both tables.  The average ES is calcualated and a CI selected by taking the min of the lower CI and the max of the upper CI from the two separate calculations.  This would be a conservative estimate. 
				
				CIoutput_Male=ci.pd(aa=allCountMatrices$male[2,1], bb=allCountMatrices$male[2,2], cc=allCountMatrices$male[1,1] ,dd=allCountMatrices$male[1,2],		method = "Nc",		alpha = 0.05, conf.level=0.95,		digits = 3,		print = FALSE,		detail.labs = FALSE )
				CIoutput_Female=ci.pd(aa=allCountMatrices$female[2,1], bb=allCountMatrices$female[2,2], cc=allCountMatrices$female[1,1] ,dd=allCountMatrices$female[1,2],		method = "Nc",		alpha = 0.05, conf.level=0.95,		digits = 3,		print = FALSE,		detail.labs = FALSE )
				CIoutput=c(CIoutput_Male[c(2, 4,5,6,7)],CIoutput_Female[c(2, 4,5,6,7)], (CIoutput_Male[5]+CIoutput_Female[5])/2, min(CIoutput_Male[6],CIoutput_Female[6]), max(CIoutput_Male[7],CIoutput_Female[7]) )   
				names(CIoutput)=c("Male_AbRateWT", "Male_AbRateKO", "Male_DiffAbRate", "Male_lower95CI", "Male_upper95CI","Fem_AbRateWT", "Fem_AbRateKO", "Fem_DiffAbRate", "Fem_lower95CI", "Fem_upper95CI" , "Stage1_AvDiffAbnRate", "Stage1_lowerCI", "Stage1_upperCI")
				
				#print(CIoutput)
				
				#stage 2 testing: Logistf regression comparing abnormality rates acrosses the sexes in knockout data only
				KOdatasetOnly=subset(result@analysedDataset, result@analysedDataset$Genotype!="WT")  #prepare knockout data for stage 2 function
				Stage2_InteractionTest=LR_KOdata(KOdatasetOnly, depVariableString=dependentVariable)
				
				#stage 2 testing: LR_KO alpha star p value
				LR_KOalphaStar=LR_KOalphaStar(y_m=allCountMatrices$male[2,2], y_f=allCountMatrices$female[2,2], n2_m=(allCountMatrices$male[2,2]+allCountMatrices$male[1,2]), n2_f=(allCountMatrices$female[2,2]+allCountMatrices$female[1,2]))
				names(LR_KOalphaStar)	=c("LR_KOalphaStar")					
				
				
				#stage 1 testing: MH_test
				Stage1_MH=MH_test(y_f=allCountMatrices$female[2,2], yx_f=(allCountMatrices$female[2,2]+allCountMatrices$female[2,1]), n2_f=(allCountMatrices$female[2,2]+allCountMatrices$female[1,2]), n1_f=(allCountMatrices$female[1,1]+allCountMatrices$female[2,1]), y_m =allCountMatrices$male[2,2], yx_m=(allCountMatrices$male[2,2]+allCountMatrices$male[2,1]), n2_m=(allCountMatrices$male[2,2]+allCountMatrices$male[1,2]), n1_m =(allCountMatrices$male[1,1]+allCountMatrices$male[2,1]), example = F)			
				names(Stage1_MH)=c("MH_alpha_star", "MH_pv", "Stage1_MH_midpv")			
				
				
				#stage 2 testing: zelen_test with one added rare event for females WT and one for males WT
				#' Females:             Males: 
				#'   
				#' |  |KO  |WT  |    |  |  |KO  |WT  |    |
				#' |--|:--:|:--:|:--:|  |--|:--:|:--:|:--:|
				#' |1 |y_f |    |yx_f|  |1 |y_m |    |yx_m|
				#' |0 |    |    |    |  |0 |    |    |    |
				#' |  |n2_f|n1_f|    |  |  |n2_m|n1_m|    |
				
				
				Stage2_zelen=zelen_test(y_f=allCountMatrices$female[2,2], yx_f=(allCountMatrices$female[2,2]+allCountMatrices$female[2,1]) + 1, n2_f=(allCountMatrices$female[2,2]+allCountMatrices$female[1,2]), n1_f=(allCountMatrices$female[1,1]+allCountMatrices$female[2,1]), y_m =allCountMatrices$male[2,2], yx_m=(allCountMatrices$male[2,2]+allCountMatrices$male[2,1]) + 1, n2_m=(allCountMatrices$male[2,2]+allCountMatrices$male[1,2]), n1_m =(allCountMatrices$male[1,1]+allCountMatrices$male[2,1]), example = F)			
				names(Stage2_zelen)=c("zelen_alpha_star", "zelen_pv", "zelen_midpv")	
				
				
				#Alternate stage2 testing: FE comparing abnormality rate acrosses the sexes in knockout data only
				#challenge.df  <-matrix(c(allCountMatrices$male[2,2],allCountMatrices$male[1,2],allCountMatrices$female[2,2],allCountMatrices$female[1,2]),
				#		nrow = 2,	dimnames = list(phenotype = c("1", "0"),	Sex = c("Male", "Female")))
				#Stage2_FE=fisher.test(challenge.df)$p.value
				#names(Stage2_FE)=c("FE_Stage2_KOonly")
				
				
				#CI for stage 2 test
				CI_stage2=ci.pd(aa=allCountMatrices$female[2,2], bb=allCountMatrices$male[2,2], cc=allCountMatrices$female[1,2] ,dd=allCountMatrices$male[1,2],	method = "Nc",alpha = 0.05, conf.level=0.95,		digits = 3,		print = FALSE,		detail.labs = FALSE )
				CI_stage2_output=c(CI_stage2[c(5,6,7)])   
				names(CI_stage2_output)=c("Stage2DiffAbRate", "Stage2lower95CI", "Stage2upper95CI" )
				
				
				CombinedResults=c(Stage1_MH, CIoutput,Stage2_InteractionTest,  CI_stage2_output, Stage2_zelen, LR_KOalphaStar)
				#print(CombinedResults)
				
				return(as.numeric(CombinedResults))
}						
