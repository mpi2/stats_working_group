# TODO:Sexual Dimorphism Classification Tag
# Author: nk3
###############################################################################

#test dataset
#setwd("T:/MM output/Feb13v2")
#df=read.csv("Eq1_merged.csv", header=TRUE, sep=",", dec=".", na.strings = "NA" )

ChangeClassification<-function(df, Phenotype_threshold=0.01){
	
	#generate an index so can walk down each row of data
	NoRows=length(df$Null_test_pvalue)
	df$Index=c(1:NoRows)
	numberTests=c(1:NoRows)
	results=c()
		
	for (bob in numberTests) {
		df_test=subset(x=df, df$Index==bob)
				
		#added a layer for no fitting of data ie NA
		if(is.na(df_test$Null_test_pvalue)==TRUE){
			ChangeClassification==NA
		}else if(df_test$Null_test_pvalue>Phenotype_threshold){
			ChangeClassification="no significant change"
		}else{
			if(df_test$Interaction_Included==FALSE){
				ChangeClassification="Significant change - both sexes equally"
			}else if(df_test$FvKO_pvalue>=0.05 && df_test$MvKO_pvalue>=0.05){
				ChangeClassification="Significant change - cannot classify"
			}else if(df_test$FvKO_pvalue<0.05 && df_test$MvKO_pvalue>=0.05){
				ChangeClassification="Significant change - females only"
			}else if(df_test$FvKO_pvalue>=0.05 && df_test$MvKO_pvalue<0.05){
				ChangeClassification="Significant change - males only"
			}else if(df_test$FvKO_par>0 && df_test$MvKO_par>0 |df_test$FvKO_par<0 && df_test$MvKO_par<0){
				if(abs(df_test$FvKO_par)>abs(df_test$MvKO_par)){
					ChangeClassification="Significant change - different size as females greater"  # change could be positive or negative but size change greater
				}else{
					ChangeClassification="Significant change - different size as males greater"
				}		
			}else{
				ChangeClassification="Significant change - different direction for the sexes"
			}										
		}
		results=c(results, ChangeClassification)	
	}
	return(results)
}


#dataTest3$ChangeClassification=ChangeClassification(dataTest3, Phenotype_threshold=0.01)
#dataTest3

