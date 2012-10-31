#----------------------------------------------------------------------------------------------------
#Function model_Formula:  set the basic fully loaded model formula which depends on how many genders levels exist.
colnames(data)
setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBAU")
data=read.csv("MBAU_data.csv", header=TRUE, sep=",", dec=".")
data2=read.csv("MBAU_data_maleonly.csv", header=TRUE, sep=",", dec=".")
model_Formula(data, "Lean.Mass")
model_Formula(data, "Fat.Mass")

#correct answer: Lean.Mass ~ Genotype + Gender + Genotype * Gender 
#correct answer: Fat.Mass ~ Genotype + Gender + Genotype * Gender 
model_Formula(data2, "Lean.Mass")
model_Formula(data2, "Fat.Mass")

#correct answer: Lean.Mass ~ Genotype 
#correct answer: Fat.Mass ~ Genotype 


#---------------------------------------------------------------------------------------------------------------------------
#Function final_genotype_model: testing the fixed effects and building final genotype model formula

setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBAU")
data=read.csv("MBAU_data.csv", header=TRUE, sep=",", dec=".")
data2=read.csv("MBAU_data_maleonly.csv", header=TRUE, sep=",", dec=".")

a=final_genotype_model(data,"Nose.To.Tail.Base.Length")
a
#correct output:  Nose.To.Tail.Base.Length ~ Genotype + Gender 

a=final_genotype_model(data2,"Nose.To.Tail.Base.Length")
a
#correct output: Nose.To.Tail.Base.Length ~ Genotype 
a=final_genotype_model(data,"Bone.Mineral.Density")
a
#correct output:Bone.Mineral.Density ~ Genotype + Gender

setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBJE study for AKG")
data=read.csv("DEXAdata.csv", header=TRUE, sep=",", dec=".")

a=final_genotype_model(data,"Nose.To.Tail.Base.Length")
#correct output:  Nose.To.Tail.Base.Length ~ Genotype + Gender 


#---------------------------------------------------------------------------------------------------------------------------
#Function final_null_model: testing the fixed effects and building final null model

setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBAU")
data=read.csv("MBAU_data.csv", header=TRUE, sep=",", dec=".")
data2=read.csv("MBAU_data_maleonly.csv", header=TRUE, sep=",", dec=".")

a=final_null_model(data,"Nose.To.Tail.Base.Length")
a
#correct output:  Nose.To.Tail.Base.Length ~ Gender 

a=final_null_model(data2,"Nose.To.Tail.Base.Length")
a
#correct output: Nose.To.Tail.Base.Length ~ 1

a=final_null_model(data,"Bone.Mineral.Density")
a
#correct output:Bone.Mineral.Density ~ Gender

setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBJE study for AKG")
data=read.csv("DEXAdata.csv", header=TRUE, sep=",", dec=".")

a=final_null_model(data,"Nose.To.Tail.Base.Length")
a
#correct output:  Nose.To.Tail.Base.Length ~ Gender 




#--------------------------------------------------------------------
#Function: tablelength.  function return the model length which is needed to extract the data correctly.
setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBAU")
data=read.csv("MBAU_data.csv", header=TRUE, sep=",", dec=".")
tablelength(data,"Nose.To.Tail.Base.Length")
#correct answer=3

tablelength(data,"Bone.Mineral.Content")
#correct answer=3


data2=read.csv("MBAU_data_maleonly.csv", header=TRUE, sep=",", dec=".")
tablelength(data2,"Nose.To.Tail.Base.Length")
#correct answer=2

setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBJE study for AKG")
data=read.csv("DEXAdata.csv", header=TRUE, sep=",", dec=".")
tablelength(data,"Nose.To.Tail.Base.Length")
#correct answer= 3

#---------------------------------------------------------------------------
#Function: table_content.  function return the table content which indicates where to query the table to grab information but also what can be grabbed.  

setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBAU")
data=read.csv("MBAU_data.csv", header=TRUE, sep=",", dec=".")
table_content(data,"Nose.To.Tail.Base.Length")
#correct answer= Gender_sig: TRUE;  WEight_sig: FALSE; Interaction-sig: FALSE

table_content(data,"Bone.Mineral.Content")
#correct answer= Gender_sig: TRUE;  WEight_sig: FALSE; Interaction-sig: FALSE


data2=read.csv("MBAU_data_maleonly.csv", header=TRUE, sep=",", dec=".")
table_content(data2,"Nose.To.Tail.Base.Length")
#correct answer= Gender_sig: FALSE;  WEight_sig: FALSE; Interaction-sig: FALSE

setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBJE study for AKG")
data=read.csv("DEXAdata.csv", header=TRUE, sep=",", dec=".")
table_content(data,"Nose.To.Tail.Base.Length")
#correct answer= Gender_sig: TRUE;  WEight_sig: FALSE; Interaction-sig: FALSE



#----------------------------------------------------------------------
#Function: finalmodel_info - Capturing the information from the final model

#output vector which includes the following: 
#1.	DepVariable
#2 Batch significance True/Flase
#3.	Variance significance True/Flase
#4.	Null test significance- p-value
#5.	genotype parameter estimate 
#6.	genotype s.e. estimate
#7.	Genotype effect p-value
#8.	gender parameter estimate 
#9.	gender s.e. estimate
#10.	Gender effect p-value
#11.	Interaction parameter estimate 
#12.	interaction s.e. estimate
#13.	interaction effect p-value
#14.	weight parameter estimate 
#15.	weight s.e. estimate
#16.	weight effect p-value
#17.	Gp1 genotype
#18.	Gp1 residuals normality test
#19.	Gp2 genotype
#20.	Gp2 residuals normality test
#21.	Blups test 
#22.	Rotated residuals normality test

setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBAU")
data=read.csv("MBAU_data.csv", header=TRUE, sep=",", dec=".")

finalmodel_info(data,"Nose.To.Tail.Base.Length")
#example when interaction not significant but weight and gender are
#correct output:
#[1] "TRUE"                 "TRUE"                 "1.08167585999297e-10"
# [4] "-0.490724060991474"   "0.072940396326446"    "5.7041801096953e-11" 
# [7] "0.372787168323284"    "0.0196571128851033"   "2.77806715040758e-58"
#[10] "NA"                   "NA"                   "NA"                  
#[13] "NA"                   "NA"                   "NA"                  
#[16] "+/+"                  "0.0617487598332429"   "Sparc/Sparc"         
#[19] "0.488295207292955"    "0.526586555557937"    "0.321734340114324"   

data2=read.csv("MBAU_data_maleonly.csv", header=TRUE, sep=",", dec=".")
finalmodel_info(data2,"Nose.To.Tail.Base.Length")
#example only one gender
#correct output
#[1] "TRUE"                 "TRUE"                 "2.66155146245728e-05"
 #[4] "-0.496505204208034"   "0.111015531821935"    "6.2494726969724e-05" 
 #[7] "NA"                   "NA"                   "NA"                  
#[10] "NA"                   "NA"                   "NA"                  
#[13] "NA"                   "NA"                   "NA"                  
#[16] "+/+"                  "0.105075344570016"    "Sparc/Sparc"         
#[19] "0.317834159944304"    "0.728186529004021"    "0.783060357900962"   
 

setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBJE study for AKG")
data=read.csv("DEXAdata.csv", header=TRUE, sep=",", dec=".")
finalmodel_info(data,"Nose.To.Tail.Base.Length")
#example all factors significant
#correct output
#"TRUE"                 "TRUE"                 "0.000708322566295827"
# [4] "-0.363244598675071"   "0.104576770965955"    "0.000862117382569142"
# [7] "0.261695181945918"    "0.0442656383892402"   "9.66960581779434e-08"
#[10] "NA"                   "NA"                   "NA"                  
#[13] "NA"                   "NA"                   "NA"                  
#[16] "+/+"                  "0.275530228756749"    "Slc25a21/Slc25a21"   
#[19] "0.524721636435985"    "0.0908468254376326"   "0.664250086033288"   

a= final_genotype_model(data, "Nose.To.Tail.Base.Length")
a

setwd("Z:/Natasha Karp files from 2008-2012/Data from tash personal drive/A New approach to filing 2010 on/Mixed model approach to analysis/sept2011/MBAU")
data=read.csv("MBAU_data.csv", header=TRUE, sep=",", dec=".")
finalmodel_info(data,"Bone.Mineral.Content")

finalmodel(data,"Bone.Mineral.Content")
#example of set where weight significant but nothing else
#correct values
#[1] "TRUE"                 "TRUE"                 "3.08946041915892e-06"
# [4] "-0.0595609242852041"  "0.0125153516148337"   "2.69663853192837e-06"
# [7] "0.0676205741356626"   "0.00416852677860946"  "4.07797629756809e-46"
#[10] "NA"                   "NA"                   "NA"                  
#[13] "NA"                   "NA"                   "NA"                  
#[16] "+/+"                  "0.588451644046573"    "Sparc/Sparc"         
#[19] "0.764222412552966"    "0.497508774683462"    "0.9830129268004"     


summary(finalmodel(data,"Bone.Mineral.Content" ))

