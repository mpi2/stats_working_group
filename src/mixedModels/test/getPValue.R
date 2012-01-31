# appends argument term to the end of formula object
addtermToFormula <- function(formulaObject, term) {
  formulaObject[[3]] <- call("+", formulaObject[[3]], substitute(term))
  return(formulaObject)
}

data=read.csv('/home/hughm/Documents/data/mixedModels/implementation/dexadataClk.csv', header=TRUE, sep=",", dec=".")
colnames(data)
attach(data)
library(nlme)

# construct batch and no batch models
loadedMeanModel=lme(Nose.To.Tail.Base.Length~Genotype + Gender + Weight + Genotype*Gender, random=~1|Assay.Date, data, na.action="na.omit", method="REML")
noBatchModel=gls(Nose.To.Tail.Base.Length~Genotype + Gender + Weight + Genotype*Gender, data, na.action="na.omit") 
# Work out if they are significantly different
canOmittBatchEffects = (anova(loadedMeanModel, noBatchModel)$p[2] / 2)> 0.05

if(!canOmittBatchEffects) {  # if batch effects are significant
  # use the batch model (linear mixed model?) for the subsequent tests
  meanModelToUse = loadedMeanModel
  # construct a similar model with unequal variance
  heterovariantModel=lme(Nose.To.Tail.Base.Length~Genotype + Gender + Weight + Genotype*Gender, random=~1|Assay.Date, data, na.action="na.omit", method="REML", weights=varIdent(form=~1|Genotype))
} else {
  # use the no batch model (Generalized Least Squares?) for the subsequent tests
  meanModelToUse = noBatchModel
  # construct a similar model with unequal variance
  heterovariantModel=gls(Nose.To.Tail.Base.Length~Genotype + Gender + Weight + Genotype*Gender, data, na.action="na.omit", weights=varIdent(form=~1|Genotype))
}

# work out if the unequal variance is significantly different
canAssumeEqualVariance = anova(meanModelToUse, heterovariantModel)$p[2] > 0.05

if(canAssumeEqualVariance) { # if variances are similar
  # use the equal variance model
  modelForReduction = meanModelToUse;
} else {
  # use the unequal variance model
  modelForReduction = heterovariantModel
}
# work out which of the variables are significant
resultOfReductionAnalysis = anova(heterovariantModel)

keepGender = resultOfReductionAnalysis$p[3] > 0.05
keepWeight = resultOfReductionAnalysis$p[4] > 0.05
keepGenotypeGenderInteraction = resultOfReductionAnalysis$p[5] > 0.05

# assemble formula, adding terms if they are significant
formulaGenotype = Nose.To.Tail.Base.Length~Genotype + Gender
formulaEffectSize = Nose.To.Tail.Base.Length~Genotype

if(keepGender) {
  formulaGenotype = addtermToFormula(formulaGenotype, Gender)
  formulaEffectSize = addtermToFormula(formulaNull, Gender)
}
if(keepWeight) {
  formulaGenotype = addtermToFormula(formulaGenotype, Weight)
  formulaEffectSize = addtermToFormula(formulaNull, Weight)
}
if(keepGenotypeGenderInteraction) {
  formulaGenotype = addtermToFormula(formulaGenotype, Genotype*Gender)
  formulaEffectSize = addtermToFormula(formulaNull, Genotype*Gender)
}

if(keepGender) {
  if(keepWeight) {
    formulaNull = Nose.To.Tail.Base.Length~Gender + Weight
  } else {
    formulaNull = Nose.To.Tail.Base.Length~Gender
  }
} else {
  if(keepWeight) {
    formulaNull = Nose.To.Tail.Base.Length~Weight
  } else {
    formulaNull = Nose.To.Tail.Base.Length~1
  }
}

# construct final models depending on whether batch effects can be omitted and if there is equal variance
if(canAssumeEqualVariance) {
  if(canOmittBatchEffects) {
    finalModelGenotype = gls(formulaGenotype , data, na.action="na.omit", method="ML")
    finalModelNull = gls(formulaNull , data, na.action="na.omit", method="ML")
    finalModelEffectSizes = gls(formulaEffectSize , data, na.action="na.omit", method="REML")
  } else {
    finalModelGenotype = lme(formulaGenotype, random=~1|Assay.Date, data, na.action="na.omit", method="ML")
    finalModelNull = lme(formulaNull, random=~1|Assay.Date, data, na.action="na.omit", method="ML")
    finalModelEffectSizes = lme(formulaEffectSize, random=~1|Assay.Date, data, na.action="na.omit", method="REML")
  }
} else {
  if(canOmittBatchEffects) {
    finalModelGenotype = gls(formulaGenotype , data, na.action="na.omit", method="ML", weights=varIdent(form=~1|Genotype))
    finalModelNull = gls(formulaNull , data, na.action="na.omit", method="ML", weights=varIdent(form=~1|Genotype))
    finalModelEffectSizes = gls(formulaEffectSize , data, na.action="na.omit", method="REML", weights=varIdent(form=~1|Genotype))
  } else {
    finalModelGenotype = lme(formulaGenotype, random=~1|Assay.Date, data, na.action="na.omit", method="ML", weights=varIdent(form=~1|Genotype))
    finalModelNull = lme(formulaNull, random=~1|Assay.Date, data, na.action="na.omit", method="ML", weights=varIdent(form=~1|Genotype))
    finalModelEffectSizes = lme(formulaEffectSize, random=~1|Assay.Date, data, na.action="na.omit", method="REML", weights=varIdent(form=~1|Genotype))
  }  
}

print('final p value = ')
print(anova(finalModelGenotype, finalModelNull)$p[2])

print(finalModelEffectSizes)
print('Something like effect size')
print(finalModelEffectSizes$coefficients[[1]][[2]])

print('canOmittBatchEffects')
print(canOmittBatchEffects)
print('canAssumeEqualVariance')
print(canAssumeEqualVariance)
print('keepGender')
print(keepGender)
print('keepWeight')
print(keepWeight)
print('keepGenotypeGenderInteraction')
print(keepGenotypeGenderInteraction)

