#test set 1
#set working directory ie where R looks, load data and attach data such the columns within the dataframe can be referenced directly
setwd("Z:/Natasha Karp files from 2008-2010/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBAU")
data=read.csv("MBAU_data.csv", header=TRUE, sep=",", dec=".")
attach(data)

#test set 2
detach(data)
setwd("Z:/Natasha Karp files from 2008-2010/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBKA")
data=read.csv("DEXAdata_homonly.csv", header=TRUE, sep=",", dec=".")
attach(data)
#het removed as too few for the analysis.

#test set 3
detach(data)
setwd("Z:/Natasha Karp files from 2008-2010/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MAAJ CC study")
data=read.csv("MAAJ_CCdata_noHET.csv", header=TRUE, sep=",", dec=".")
attach(data)

#test set 4 - HOM's and HET's and WT
setwd("Z:/Natasha Karp files from 2008-2010/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBKA")
data=read.csv("DEXAdata.csv", header=TRUE, sep=",", dec=".")
attach(data)

library(nlme)

#Note:
#assumes data has been cleaned to remove genotype failed such that you only have WT, homs or hets.  
#Haven't considered Hemi classification - Need to know how we want hemi analysed - treated as hom's?  If so they need recoding to hom's. 
#what have I not wrapped? The input of the data and the loading of the library(nlme) and library(car),  attachment and detachment of data. 
#naming issues will be an issue.  For example,  body weight is sometimes denoted by Weight sometimes Weight.Value.  Currently I have based on the reporting headers but the database uses different headers.

#What about lines which have both: Hom's and Het's?
#You can analyse them simultaneously,  however if the n is too low for a genotype group eg one hom there is no value to including (eg MAAJ).  I have currently been manually curating the data for this. I have not written a function for this.  
#Does the code work if analysed simultaneously? Does the output make sense?  - Yes

#-----------------------------------------------------
#Function1: function to test batch   
# A full loaded model formula (ie one with all variables of interest) is assembled.  Weight cannot be in the model as a dependent variable if it is the independent variable - hence the if else questions to build the model. 
# Then the model formula are used to build two models one with batch included as a random effect and hence uses a mixed model and one where batch is not included and uses linear regression
# The two models are compared via anova to test the null hypothesis that batch is not significant (Ho) where the alternative hypothesis (Ha) is that batch is significant.  If the p value <0.05,  then the we reject the null hypothesis and accept the alternative and say we need to keep batch  
# the p-value for batch test is divided by 2 as the test is on the boundary of the parameter space for a variance and the null distribution of the likelihood ratio test statistics follows a mixture of chi-squared distributions with equal weight of 0.5.  
# As the output is yes or no on keep_batch then if p<0.05,  we say variance should be equal

testing_batch<-function(dataset, depVariable){
    if(depVariable=="Weight"){
      model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Genotype*Gender", sep= "+")))
      }
    else  {
      model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Weight", "Genotype*Gender", sep= "+")))
      }
    model_withbatch=do.call("lme", args = list(model.formula, random=~1|Assay.Date, dataset, na.action="na.omit", method="REML"))
    model_withoutbatch=do.call("gls", args=list(model.formula, dataset, na.action="na.omit"))
    p.value=(anova(model_withbatch, model_withoutbatch)$p[2])/2
    keep_batch= p.value<0.05
    return(keep_batch)
 
 }
#test of batch function
a=testing_batch(data, "Nose.To.Tail.Base.Length")
a
a=testing_batch(data, "Lean.Mass")
a
a=testing_batch(data, "Weight")
a
a=testing_batch(data, "Tp")
a
a=testing_batch(data, "Na")
a
a=testing_batch(data, "Glyc")
a
#---------------------------------------------------
#Function 2: function to test equal variance

# A full loaded model forumula(ie one with all variables of interest) is assembled.  Weight cannot be in the model as a dependent variable if it is the independent variable - hence the if else questions to build the model. 
# Then the models formula are used to build two models one with homogenous variance and one with heterogenous variance
# The two models are compared via anova to test the null hypothesis that variance is equal (Ho) where the alternative hypothesis (Ha) is that variance is not equal.  If the p value <0,05,  the null hypotehsis is rejected and the alternative accepted. 
# As the output is yes or no on keep_equalvar then if p<0.05,  we say variance should be equal

testing_variance<-function(dataset, depVariable){
    
    if(depVariable=="Weight"){
      model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Genotype*Gender", sep= "+")))
    }else {
      model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Weight", "Genotype*Gender", sep= "+")))
    }
   model_equalvariance=do.call("lme", args=list(model.formula, random=~1|Assay.Date, dataset, na.action="na.omit", method="REML") )
   model_hetvariance= do.call("lme", args=list(model.formula, random=~1|Assay.Date, dataset, weights=varIdent(form=~1|Genotype), na.action="na.omit", method="REML"))
   p.value=(anova(model_equalvariance, model_hetvariance)$p[2])
   keep_equalvar= p.value>0.05
   return(keep_equalvar)
}
#test of variance function
a=testing_variance(data, depVariable="Nose.To.Tail.Base.Length")
a
a=testing_variance(data, depVariable="Lean.Mass")
a
a=testing_variance(data, depVariable="Weight")
a
a=testing_variance(data, depVariable="Tp")
a
a=testing_variance(data, depVariable="Na")
a
a=testing_variance(data, depVariable="Glyc")
a

#------------------------------------------------------------------
# Function 3: function to build model following variance and batch test for testing of fixed effects
# Goal: To use the results of F1 (testing_batch) and F2 (testing_variance) a model is built
# Weight cannot be in the model as a dependent variable if it is the independent variable - hence the if else questions to build the model. 
# The keep_batch and keep_equal variance functions are called to build the correct models via a series of if else rules. 

model_forFIXEDtest<-function(dataset, depVariable){

    #test whether batch is significant
    keep_batch=testing_batch(dataset,depVariable)

    #test whether variance is significant
    keep_equalvar=testing_variance(dataset, depVariable)

    #build correct model depending on the output of the two questions (batch and variance)
   
    if(depVariable=="Weight"){
      model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Genotype*Gender", sep= "+")))
    }else{    
      model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Weight", "Genotype*Gender", sep= "+")))
    }  
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


#test of function bulding the correct model
t=model_forFIXEDtest(data,  "Nose.To.Tail.Base.Length")
t
t=model_forFIXEDtest(data,  "Lean.Mass")
t
t=model_forFIXEDtest(data,  "Weight")
t
t=model_forFIXEDtest(data,  "Tp")
t
t=model_forFIXEDtest(data,  "Na")
t
t=model_forFIXEDtest(data,  "Glyc")
t
#---------------------------------------------------------------------------------------------------------------------------
#Function 4: testing the fixed effects and building final genotype model formula
# Goal:  to test fixed effects of the model and based on the output build the final genotype model formula for later testing.  As a genotype model it automatically includes genotype. 
# First F3 is called (model_forFIXEDtest) to build the model to be queried 
# The anova function tests the fixed effects associated by treatment with a null hypothesis that the regression coefficients are equal to zero  and an alternative hypothesis that the regression coefficient are not equal to zero.
# If the p-values of these tests are less than 0.05 we reject the null and accept the alternative that the are significant components of the model and should be included. 
# Note a complexity surrounds the interaction term  - if it is significant but gender is excluded it is excluded. 
# Weight cannot be in the model as a dependent variable if it is the independent variable - hence the if else questions to build the formula. 
  

final_genotype_model<-function(dataset, depVariable){
  model_followingbatchandvartest=model_forFIXEDtest(dataset, depVariable)
  anova_results = anova(model_followingbatchandvartest, type="marginal")$"p-value" < 0.05
  keepGender = anova_results[3]
  
   
    if(depVariable=="Weight"){
      #test whether fixed effects of the model are significant
      keepInteraction = anova_results[4]
      if(keepGender && keepInteraction){
        return(model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Genotype*Gender", sep= "+"))))
      }else if(keepGender &&  !keepInteraction){
        return(model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender", sep= "+"))))
      }else if(!keepGender &&  !keepInteraction){
        return(model.formula <- as.formula(paste(depVariable, "~", "Genotype")))
      }else if(!keepGender  && keepInteraction){
        return(model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Genotype*Gender", sep= "+"))))
      } 
    }else{
      keepWeight = anova_results[4]
      keepInteraction = anova_results[5]
      if ((keepGender && keepWeight && keepInteraction)| (!keepGender && keepWeight && keepInteraction)){
          return(model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Weight", "Genotype*Gender", sep= "+"))))
      } else if((keepGender && !keepWeight && keepInteraction)|(!keepGender && !keepWeight && keepInteraction)){
        return(model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender",  "Genotype*Gender", sep= "+"))))
      } else if(keepGender && !keepWeight && !keepInteraction){
        return(model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender", sep= "+"))))
      } else if(!keepGender && !keepWeight && !keepInteraction){
        return(model.formula <- as.formula(paste(depVariable, "~", "Genotype")))
      } else if(keepGender && keepWeight && !keepInteraction){
        return(model.formula <- as.formula(paste(depVariable, "~", paste("Genotype", "Gender", "Weight", sep= "+"))))
      } else if(!keepGender && keepWeight && !keepInteraction){
        return(model.formula <- as.formula(paste(depVariable, "~",paste("Genotype", "Weight", sep= "+"))))
      }
   }
}

b=final_genotype_model(data,"Nose.To.Tail.Base.Length")
b
b=final_genotype_model(data,"Lean.Mass")
b
b=final_genotype_model(data,"Weight")
b
b=final_genotype_model(data,"Tp")
b
b=final_genotype_model(data,"Na")
b
b=final_genotype_model(data,"Glyc")
b
#---------------------------------------------------------------------------------------------------------------------------
#Function 5: testing the fixed effects and building final null model
# Goal:  to test fixed effects of the model and based on the output build the final null model formula for later testing - as a null model it automatically excludes genotype and interaction term. 
# First F3 is called (model_forFIXEDtest) to build the model to be queried 
# The anova function tests the fixed effects associated by treatment with a null hypothesis that the regression coefficients are equal to zero  and an alternative hypothesis that the regression coefficient are not equal to zero.
# If the p-values of these tests are less than 0.05 we reject the null and accept the alternative that the are significant components of the model and should be included. 
# Weight cannot be in the model as a dependent variable if it is the independent variable - hence the if else questions to build the formula. 
# If no terms are significant a model can be build with just an intercept element this is specified as  "model.formula <- as.formula(paste(depVariable, "~", "1"))"


final_null_model<-function(dataset, depVariable){ 
    model_followingbatchandvartest=model_forFIXEDtest(dataset, depVariable)
    anova_results = anova(model_followingbatchandvartest, type="marginal")$"p-value" < 0.05
    keepGender = anova_results[3]
    
    if(depVariable=="Weight"){
      #test whether fixed effects of the model are significant
     
      keepInteraction = anova_results[4]
              
      if(!keepGender && !keepInteraction){
        return(model.formula <- as.formula(paste(depVariable, "~", "1")))
      }else{
        return(model.formula <- as.formula(paste(depVariable, "~", "Gender")))
      } 
    }else{
        
        keepWeight = anova_results[4]
        keepInteraction = anova_results[5]
        
       if(!keepGender && !keepWeight && !keepInteraction){
        return(model.formula <- as.formula(paste(depVariable, "~", "1")))
       }else if (!keepGender && keepWeight && !keepInteraction){
        return(model.formula <- as.formula(paste(depVariable, "~", "Weight")))
       }else if((keepGender && !keepWeight && !keepInteraction)|(!keepGender && !keepWeight && keepInteraction)|(keepGender && !keepWeight && keepInteraction)){
        return(model.formula <- as.formula(paste(depVariable, "~", "Gender")))
       }  
       else{
        return(model.formula <- as.formula(paste(depVariable, "~", paste("Gender", "Weight", sep= "+"))))
       }       
    }
}

b=final_null_model(data,"Nose.To.Tail.Base.Length")
b
b=final_null_model(data,"Lean.Mass")
b
b=final_null_model(data,"Weight")
b
b=final_null_model(data,"Tp")
b
b=final_null_model(data,"Na")
b
b=final_null_model(data,"Glyc")
b

#----------------------------------------------------------------------------------
#Function 6: testing the genotype effect
# The goal of the function is to give a pvalue of whether genotype effect is significant by compare the genotype model with the null model with the anova function
# The model_formula_null and model_formula_genotype are called to define the models for comparison. 
# The keep_batch and keep_equal variance function are called this allows if but rules to build the correct model comparison ie if batch is included then we use a mixed model etc
# For each possible combination,  then an anova model is used to report the pvalue.
# For testing the genotype effect we use method= ML


testing_genotype_effect<-function(dataset, depVariable){
    
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
        model_genotype=do.call("gls", args = list(model_formula_genotype,  dataset,weights=varIdent(form=~1|Genotype), na.action="na.omit"))
        model_null=do.call("gls", args=list(model_formula_null, dataset,weights=varIdent(form=~1|Genotype), na.action="na.omit"))
        return(p.value=(anova(model_genotype, model_null)$p[2]))
    }else if(!keep_batch && keep_equalvar){
        model_genotype=do.call("gls", args = list(model_formula_genotype,  dataset, na.action="na.omit"))
        model_null=do.call("gls", args=list(model_formula_null, dataset, na.action="na.omit"))
        return(p.value=(anova(model_genotype, model_null)$p[2]))
    } 
}  
    
d=testing_genotype_effect(data,"Nose.To.Tail.Base.Length")
d
d=testing_genotype_effect(data,"Lean.Mass")
d
d=testing_genotype_effect(data,"Weight")
d
d=testing_genotype_effect(data,"Tp")
d
d=testing_genotype_effect(data,"Na")
d
d=testing_genotype_effect(data,"Glyc")
d
#-----------------------------------------------------------------------------------------
#Function 7: following function returns the final model which is needed for the diagnostic plots 
#Goal of this function is to return a model output that is the final model following all the previous tests which can be used to generate diagnostics plots and output the final model details. 
# The model_formula_null and model_formula_genotype are called to define the models for comparison. 
# The keep_batch and keep_equal variance function are called this allows if but rules to build the correct model comparison ie if batch is included then we use a mixed model etc
# For each possible combination,  then the model is fitted to the data and the model reported as the output.
#this is a separate function to above as for the estimates of the fixed effects we use method=REML whilst for the genotype test we used method= ML   
#na.action = na.exclude as in this form the residue calculated from the model output have the same length as the original datafile
#see http://stackoverflow.com/questions/6882709/how-do-i-deal-with-nas-in-residuals-in-a-regression-in-r
      

finalmodel<-function(dataset, depVariable){

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

finalmodel(data,"Nose.To.Tail.Base.Length")
finalmodel(data,"Lean.Mass")
finalmodel(data,"Weight")
finalmodel(data,"Tp")
finalmodel(data,"Na")
finalmodel(data,"Glyc")
 
#--------------------------------------------------------------------
#Function 8:  following function return the model length which is needed to extract the data correctly. 
#Goal of this function is to return a length measure that will determine the table width so when we try to capture the output we can point correctly to the right point in a vector of values. 


 tablelength<-function(dataset,depVariable){
   model_followingbatchandvartest=model_forFIXEDtest(dataset, depVariable)
   anova_results = anova(model_followingbatchandvartest, type="marginal")$"p-value" < 0.05
   keepGender = anova_results[3]
   
       
   #here I am counting how many elements will be on the table based on how the model used.  However the model used isnot a direct additive system of the yes calls (ie you can't have an interaction if you don't have gender in the model).    
      if(depVariable=="Weight"){
         keepInteraction = anova_results[4]
         if((keepGender && keepInteraction)|(!keepGender && keepInteraction)){
         return(table_length=4)
         }else if(keepGender && !keepInteraction){
         return(table_length=3)
         }else {
         return(table_length=2)
               }
     
      }else{
         keepWeight = anova_results[4]
         keepInteraction = anova_results[5]
           if((keepGender && keepInteraction && keepWeight)|(!keepGender && keepInteraction && keepWeight)){
           return(table_length=5)
           }else if((keepGender && keepInteraction && !keepWeight)|(keepGender && !keepInteraction && keepWeight)|(!keepGender && keepInteraction && !keepWeight)){
           return(table_length=4)
           }else if(!keepGender && !keepInteraction && !keepWeight){
           return(table_length=2)
           }else {
           return(table_length=3)
           }        
      }        
 }


summary(finalmodel(data,"Nose.To.Tail.Base.Length"))

b=final_null_model(data,"Nose.To.Tail.Base.Length")
b
e=tablelength(data,"Nose.To.Tail.Base.Length")
e
b=final_null_model(data,"Lean.Mass")
b
e=tablelength(data,"Lean.Mass")
e
b=final_null_model(data,"Weight")
b
e=tablelength(data,"Weight")
e
e=tablelength(data,"Tp")
e
e=tablelength(data,"Na")
e

e=tablelength(data,"Glyc")
e

#----------------------------------------------------------------------
#Function 9: Capturing the information from the final model
#problem with length correction

finalmodel_info<-function(dataset, depVariable){
    
    modeloutput=summary(finalmodel(dataset, depVariable))
    keep_batch=testing_batch(dataset,depVariable)
    Nulltest_genotype_pvalue=testing_genotype_effect(dataset,depVariable)
    lengthoftable=tablelength(dataset,depVariable)
    numberofGenotypes=nlevels(Genotype)
    a=levels(Genotype)

   #problem depending on the length of the table where we grab values   
       
    if(keep_batch && numberofGenotypes==2){
     #for mixed model and with two genotype comparison
        
        genotypeestimate = modeloutput[["tTable"]][[2]]
        genotypeestimate_SE = modeloutput[["tTable"]][[(2+lengthoftable)]]
        degreesoffreedom =  modeloutput[["tTable"]][[(2+2*lengthoftable)]]
        genotype_t_value =  modeloutput[["tTable"]][[(2+3*lengthoftable)]]
        genotype_p_value =  modeloutput[["tTable"]][[(2+4*lengthoftable)]]
        valuestocollect = c("genotype", "effect", "std error", "DF", "t-value", "p-value","null p_value")  
        genotype_values=c(a[2], genotypeestimate,genotypeestimate_SE,  degreesoffreedom,genotype_t_value,  genotype_p_value,Nulltest_genotype_pvalue)
        values=data.frame(valuestocollect, genotype_values)
        return(values)
        
    }else if(keep_batch && numberofGenotypes==3){
      #for mixed model and with three genotype comparison   add one to the length of the table
        lengthoftable=lengthoftable+1
        genotypeestimate = modeloutput[["tTable"]][[2]]
        genotypeestimate_SE = modeloutput[["tTable"]][[(2+lengthoftable)]]
        degreesoffreedom =  modeloutput[["tTable"]][[(2+2*lengthoftable)]]
        genotype_t_value =  modeloutput[["tTable"]][[(2+3*lengthoftable)]]
        genotype_p_value =  modeloutput[["tTable"]][[(2+4*lengthoftable)]]
        valuestocollect = c("genotype","effect", "std error", "DF", "t-value", "p-value","null p_value")  
        genotype_values=c(a[2], genotypeestimate,genotypeestimate_SE,  degreesoffreedom,genotype_t_value,  genotype_p_value,Nulltest_genotype_pvalue)
        values=data.frame(valuestocollect, genotype_values)
         genotypeestimate2 = modeloutput[["tTable"]][[3]]
        genotypeestimate_SE2 = modeloutput[["tTable"]][[(3+lengthoftable)]]
        degreesoffreedom2 = modeloutput[["tTable"]][[(3+2*lengthoftable)]]
        genotype_t_value2 = modeloutput[["tTable"]][[(3+3*lengthoftable)]]
        genotype_p_value2 = modeloutput[["tTable"]][[(3+4*lengthoftable)]]
       genotype_values2=c(a[3], genotypeestimate2,genotypeestimate_SE2,  degreesoffreedom2,genotype_t_value2,  genotype_p_value2,Nulltest_genotype_pvalue)
        values=data.frame(valuestocollect, genotype_values, genotype_values2)  
        return(values)
    
    }else if (!keep_batch && numberofGenotypes==2){  
        #adaption for being a linear model rather than a mixed model with two genotype comparison
        genotypeestimate = modeloutput[["tTable"]][[2]]
        genotypeestimate_SE = modeloutput[["tTable"]][[(2+lengthoftable)]]
        degreesoffreedom =  "NA"
        genotype_t_value =  modeloutput[["tTable"]][[(2+2*lengthoftable)]]
        genotype_p_value =  modeloutput[["tTable"]][[(2+3*lengthoftable)]]
        valuestocollect = c("genotype", "effect", "std error", "DF", "t-value", "p-value","null p_value")  
        genotype_values=c(a[2], genotypeestimate,genotypeestimate_SE,  degreesoffreedom,genotype_t_value,  genotype_p_value,Nulltest_genotype_pvalue)
        values=data.frame(valuestocollect, genotype_values)
        return(values)
    
    }else if (!keep_batch && numberofGQenotypes==3){
       #capture for linear model with a three genotype comparison
        lengthoftable=lengthoftable+1
        genotypeestimate = modeloutput[["tTable"]][[2]]
        genotypeestimate_SE = modeloutput[["tTable"]][[(2+lengthoftable)]]
        degreesoffreedom =  "NA"
        genotype_t_value =  modeloutput[["tTable"]][[(2+2*lengthoftable)]]
        genotype_p_value =  modeloutput[["tTable"]][[(2+3*lengthoftable)]]
        valuestocollect = c("genotype","effect", "std error", "DF", "t-value", "p-value","null p_value")  
        genotype_values1=c(a[2], genotypeestimate,genotypeestimate_SE,  degreesoffreedom,genotype_t_value,  genotype_p_value,Nulltest_genotype_pvalue)
        genotypeestimate2 = modeloutput[["tTable"]][[3]]
        genotypeestimate_SE2 = modeloutput[["tTable"]][[(3+lengthoftable)]]
        degreesoffreedom2 =  "NA"
        genotype_t_value2 =  modeloutput[["tTable"]][[(3+2*lengthoftable)]]
        genotype_p_value2 =  modeloutput[["tTable"]][[(3+3*lengthoftable)]]
       genotype_values2=c(a[3], genotypeestimate2,genotypeestimate_SE2,  degreesoffreedom2,genotype_t_value2,  genotype_p_value2,Nulltest_genotype_pvalue)
        values=data.frame(valuestocollect, genotype_values, genotype_values2)
        return(values)  
     }
}

values=finalmodel_info(data,"Nose.To.Tail.Base.Length")          
values    
values=finalmodel_info(data,"Lean.Mass")
values    
values=finalmodel_info(data,"Weight")
values    
values=finalmodel_info(data,"Tp")
values   
values=finalmodel_info(data,"Na")
values   
values=finalmodel_info(data,"Glyc")
values   
#-------------------------------------------------------------------------------------------
#If data is going to be attached via wrapper function then the diagnostic plots will need the data being attached for them within their code.  depends on implementation detail. 

#for  diagnostic plots
#Function D1:  Raw data boxplot split by gender and genotype       
boxplot_bygender<-function(dataset,depVariable,v_name){
  Male = subset(dataset, Gender=="Male")
  Female= subset(dataset, Gender=="Female")
  split.screen(c(1,2), erase = TRUE)
  screen(n = 1 , new = TRUE)
  boxplot(Male[ , depVariable]~Male$Genotype, ylab=v_name, xlab="Genotype")
  legend("topright", "Male", cex=1.3, bty="n")
  close.screen(n=1, all.screens = FALSE)
  screen(n = 2 , new = TRUE)
  boxplot(Female[ , depVariable]~Female$Genotype, ylab=v_name, xlab="Genotype")
  legend("topright", "Female", cex=1.3, bty="n")
  close.screen(all.screens = TRUE)
}
boxplot_bygender(data,  "Nose.To.Tail.Base.Length", "Body length")

#---------------------------------------------------------------------------
#Function D2: QQ Residue plots for each genotype for each gender       
ResidueQQplot_bygenotype<-function(dataset, depVariable){
    modeloutput=finalmodel(dataset,depVariable)
    res=resid(modeloutput)
    data_all= data.frame(dataset, res)
    a=levels(data_all$Genotype)
    genotype_no=length(a)
    
    if(genotype_no==2){
      split.screen(c(1,2), erase = TRUE)
      screen(n = 1 , new = TRUE)
      Gp1 = subset(data_all, Genotype==a[1])
      Gp2 = subset(data_all, Genotype==a[2])
      qqnorm(Gp1$res, main=" ")
      qqline(Gp1$res)
      legend("topleft", a[1], cex=1.3, bty="n")
      close.screen(n=1, all.screens = FALSE)
      qqnorm(Gp2$res, main=" ")
      qqline(Gp2$res)
      legend("topleft", a[2], cex=1.3, bty="n")
      screen(n = 2 , new = TRUE)
      close.screen(all.screens = TRUE)
    }else{
      split.screen(c(1,3), erase = TRUE)
      screen(n = 1 , new = TRUE)
      Gp1 = subset(data_all, Genotype==a[1])
      Gp2 = subset(data_all, Genotype==a[2])
      Gp3 = subset(data_all,Genotype==a[3])
      qqnorm(Gp1$res, main=" ")
      qqline(Gp1$res)
      legend("topleft", a[1], cex=1.3, bty="n")
      close.screen(n=1, all.screens = FALSE)
      screen(n = 2 , new = TRUE)
      qqnorm(Gp2$res, main=" ")
      qqline(Gp2$res)
      legend("topleft", a[3], cex=1.3, bty="n")
      close.screen(n=2, all.screens = FALSE)
      screen(n = 3 , new = TRUE)
      qqnorm(Gp3$res, main=" ")
      qqline(Gp3$res)
      legend("topleft", a[3], cex=1.3, bty="n")
      close.screen(all.screens = TRUE)
    }
}
ResidueQQplot_bygenotype(data,  "Nose.To.Tail.Base.Length")
ResidueQQplot_bygenotype(data,  "Lean.Mass")
ResidueQQplot_bygenotype(data,  "Weight")
#-------------------------------------------------------------------------
#http://stackoverflow.com/questions/6882709/how-do-i-deal-with-nas-in-residuals-in-a-regression-in-r
#Function D3: Predicted versus residue plots split by genotype   
ResVpred_bygenotype<-function(dataset,depVariable){
    modeloutput=finalmodel(dataset,depVariable)
    pred = predict(modeloutput)
    res=resid(modeloutput)
     data_All= data.frame(dataset, res, pred)
       
   
    a=levels(dataset$Genotype)
    genotype_no=length(a)
    
    if(genotype_no==2){
      split.screen(c(1,2), erase = TRUE)
      screen(n = 1 , new = TRUE)
      Gp1pred = subset(data_All, Genotype==a[1])
      Gp2pred = subset(data_All, Genotype==a[2])
      plot(x=Gp1pred$pred, y=Gp1pred$res, xlab="Predicted", ylab="Residual")
      legend("topleft", a[1], cex=1.3, bty="n")
      close.screen(n=1, all.screens = FALSE)
      screen(n = 2 , new = TRUE)
      plot(x=Gp2pred$pred, y=Gp2pred$res, xlab="Predicted", ylab="Residual")
      legend("topleft", a[2], cex=1.3, bty="n")
      close.screen(all.screens = TRUE)
    }else{
      split.screen(c(1,3), erase = TRUE)
      screen(n = 1 , new = TRUE)
      Gp1pred = subset(data_All, Genotype==a[1])
      Gp2pred = subset(data_All, Genotype==a[2])
      Gp3pred=subset(data_All, Genotype==a[3])
      plot(x=Gp1pred$pred, y=Gp1pred$res, xlab="Predicted", ylab="Residual")
      legend("topleft", a[1], cex=1.3, bty="n")
      close.screen(n=1, all.screens = FALSE)
      screen(n = 2 , new = TRUE)
      plot(x=Gp2pred$pred, y=Gp2pred$res, xlab="Predicted", ylab="Residual")
      legend("topleft", a[2], cex=1.3, bty="n")
      close.screen(n=2, all.screens = FALSE)
      screen(n = 3 , new = TRUE)
      plot(x=Gp3pred$pred, y=Gp3pred$res, xlab="Predicted", ylab="Residual")
      legend("topleft", a[3], cex=1.3, bty="n")
      close.screen(all.screens = TRUE)
     }
}


ResVpred_bygenotype(data,  "Nose.To.Tail.Base.Length")
ResVpred_bygenotype(data,  "Weight")
ResVpred_bygenotype(data,  "Lean.Mass")

debug(ResVpred_bygenotype)

#-----------------------------------------------------------------------------------
#Function D4:  Body weight versus dependent variable scatter plot

library(car)
#did not work when library call was within function.

weight_versus_depVariable<-function(dataset, depVariable){
    model.formula <- as.formula(paste(depVariable, "~", paste("Weight", "Genotype", sep= "|")))
    scatterplot(model.formula)
}
weight_versus_depVariable(data, "Nose.To.Tail.Base.Length")
weight_versus_depVariable(data, "Lean.Mass")


----------------------------------------------------------------------------------------------------
#Function D5:    QQ plot of Blups 

blups_plot<-function(dataset, depVariable){ 
    modeloutput=finalmodel(dataset,depVariable)
    blups=ranef(modeloutput)
    qqnorm(blups[ ,1])
    qqline(blups[ ,1])
}

blups_plot(data, "Nose.To.Tail.Base.Length")
blups_plot(data, "Lean.Mass")
blups_plot(data, "Weight")

----------------------------------------------------------------------------------------------------------
#Function D6: Residue versus batch split by genotype

ResidueVersusBatch_bygenotype<-function(dataset, depVariable){
    modeloutput=finalmodel(dataset,depVariable)
    res=resid(modeloutput)
    data_all= data.frame(dataset, res)
    a=levels(dataset$Genotype)
    genotype_no=length(a)
    
    if(genotype_no==2){
    split.screen(c(1,2), erase = TRUE)
    screen(n = 1 , new = TRUE)
    Gp1 <- subset(data_all, Genotype==a[1])
    Gp2 <- subset(data_all, Genotype==a[2])
    with(Gp1, boxplot(res~Assay.Date, ylab="residues", xlab="Batch"))
    legend("bottomleft", a[1], cex=1.3, bty="n")
    screen(n = 2 , new = TRUE)
    with(Gp2, boxplot(res~Assay.Date, ylab="residues", xlab="Batch"))
    legend("bottomleft", a[2], cex=1.3, bty="n")
    close.screen(all.screens = TRUE)
     
    }else{
      split.screen(c(1,3), erase = TRUE)
      screen(n = 1 , new = TRUE)
      Gp1 <- subset(data_all, Genotype==a[1])
      Gp2 <- subset(data_all, Genotype==a[2])
      Gp3=subset(data_All, Genotype==a[3])
      with(Gp1, boxplot(res~Assay.Date, ylab="residues", xlab="Batch"))
      legend("bottomleft", a[1], cex=1.3, bty="n")
      close.screen(n=1, all.screens = FALSE)
      screen(n = 2 , new = TRUE)
      with(Gp2, boxplot(res~Assay.Date, ylab="residues", xlab="Batch"))
      legend("bottomleft", a[2], cex=1.3, bty="n")
      close.screen(n=2, all.screens = FALSE)
      screen(n = 3 , new = TRUE)
      with(Gp3, boxplot(res~Assay.Date, ylab="residues", xlab="Batch"))
      legend("bottomleft", a[3], cex=1.3, bty="n")
      close.screen(all.screens = TRUE)
             
           }
 }
    

ResidueVersusBatch_bygenotype(data, "Nose.To.Tail.Base.Length")
ResidueVersusBatch_bygenotype(data, "Lean.Mass")
ResidueVersusBatch_bygenotype(data, "Weight")

