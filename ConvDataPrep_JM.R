library(tidyverse)
library(mice)

## Import, clean, and restructure data into analysis-ready formats

### Import data prepared in feature selection process in RemFeatSelect_JM.R script
jm_long_analysis_conv = read.csv("/Users/maw/Documents/Documents/Yale/Cannon Lab/Joint Modeling Risk Pred/Conversion/Data/jm_long_analysis_11.10.22.csv")
jm_long_analysis_conv$X = NULL

jm_long_analysis_conv$VisitLabel[jm_long_analysis_conv$VisitLabel == "BL"] = 0 
jm_long_analysis_conv$VisitLabel[jm_long_analysis_conv$VisitLabel == "M2"] = 2 
names(jm_long_analysis_conv)[names(jm_long_analysis_conv) == "VisitLabel"] = "month"

jm_long_analysis_conv$months_to_conversion = jm_long_analysis_conv$days_from_baseline_until_convers/30.42
jm_long_analysis_conv$months_to_conversion = ifelse(is.na(jm_long_analysis_conv$months_to_conversion), 24, jm_long_analysis_conv$months_to_conversion)

### Import BL risk calc data 
jm_bl_analysis_conv = read.csv("/Users/maw/Documents/Documents/Yale/Cannon Lab/Risk Calc/NAPLS3/Data/NAPLS3_riskcalc_9.24.21.csv")
jm_bl_analysis_conv$X = NULL
jm_bl_analysis_conv = jm_bl_analysis_conv[,-c(6)]


## add BL variables to larger dataset - be sure to exclude variables already included in long data 
# overlap variables: P2_SOPS, P5_SOPS, CDS4
bl_variable = jm_bl_analysis_conv[,c(1,4:11)]

full_long = merge(jm_long_analysis_conv, bl_variable, by = "NAPLS_ID")

# re-order variables so they are more intuitive
full_long = full_long[,c(1,2,9,4:8,10:17)]

jm_full_conv = full_long


##### MICE #####
# MICE for missing data - only on TDP variables and baseline calc variables
jm_full_conv[,c(5:16)] = sapply(jm_full_conv[,c(5:16)], as.numeric)
jm_full_conv$converters = as.factor(jm_full_conv$converters)

## first see if there are patterns of missingness 
jm_full_conv_wide = pivot_wider(jm_full_conv,
                               names_from = "month", 
                               values_from = c(5:8))

p_missing = unlist(lapply(jm_full_conv_wide, function(x) sum(is.na(x))))/nrow(jm_full_conv_wide)
sort(p_missing[p_missing > 0], decreasing = TRUE)
# all below 25% missing


## Multiple imputation ## 
###### Multiple imputation with MICE package ######
set.seed(42)
micelong0 = mice(jm_full_conv, maxit = 0)
meth_micelong = micelong0$method
pred_micelong = micelong0$predictorMatrix

jm_imp_conv = mice(data = jm_full_conv, 
                  method = "2l.lmer",
                  pred = pred_micelong,
                  maxit = 20)

View(complete(jm_imp_conv))
jm_imp_conv_long = complete(jm_imp_conv) # computing means of all imputed data sets from MICE, storing this as new dataset

# wide format - can use this for JM analysis because it has BL variables and long variables 
jm_imp_conv_wide = pivot_wider(jm_imp_conv_long, 
                              names_from = "month",
                              names_sep = ".",
                              values_from = c(5:8))

# can use this also because the labels don't have ".0" and ".2"
jm_imp_conv_bl = filter(jm_imp_conv_long, month == 0)


### store these datasets for future use and importing 
write.csv(jm_imp_conv_long, "/Users/maw/Documents/Documents/Yale/Cannon Lab/Joint Modeling Risk Pred/Conversion/Data/JM_longdata_conversion_imputed.csv")
write.csv(jm_imp_conv_wide, "/Users/maw/Documents/Documents/Yale/Cannon Lab/Joint Modeling Risk Pred/Conversion/Data/JM_widedata_conversion_imputed.csv")
write.csv(jm_imp_conv_bl, "/Users/maw/Documents/Documents/Yale/Cannon Lab/Joint Modeling Risk Pred/Conversion/Data/JM_bldata_conversion_imputed.csv")





