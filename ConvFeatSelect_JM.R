# Import data - variables already combined, includes ~100 from clinical, demographic, neurocognitive 
mlm_data_long = read.csv("/Users/maw/Documents/Documents/Yale/Cannon Lab/Joint Modeling Risk Pred/Conversion/Data/mlm_data_long.csv")

# ensure outcome variable is a factor
mlm_data_long$converters = factor(mlm_data_long$converters)

# drop subjects who do not have SOPS data at both time points 
mlm_data_long = drop_na(mlm_data_long, cSOPSPositive)

month0 = filter(mlm_data_long, VisitLabel == "BL")
month2 = filter(mlm_data_long, VisitLabel == "M2")

month0_match = filter(month0, month0$NAPLS_ID %in% month2$NAPLS_ID)
month2_match = filter(month2, month2$NAPLS_ID %in% month0$NAPLS_ID)

mlm_long_analysis = data.frame(rbind(month0_match, month2_match))
mlm_bl_analysis = filter(mlm_long_analysis, VisitLabel == "BL")


# Remove variables (summary scores) that are colinear with individual items
mat = round(cor(mlm_bl_analysis[,c(9,17:116)], use = "pairwise.complete.obs", method = "pearson"), 2)
cor_df = mlm_bl_analysis[,c(9,17:116)]

corr_simple <- function(data = cor_df, sig = 0.4){
  #convert data to numeric in order to run correlations
  #convert to factor first to keep the integrity of the data - each value will become a number rather than turn into NA
  df_cor <- data %>% mutate_if(is.character, as.factor)
  df_cor <- df_cor %>% mutate_if(is.factor, as.numeric)
  #run a correlation and drop the insignificant ones
  corr <- cor(df_cor)
  #prepare to drop duplicates and correlations of 1     
  corr[lower.tri(corr, diag = TRUE)] <- NA 
  #drop perfect correlations
  corr[corr == 1] <- NA 
  #turn into a 3-column table
  corr <- as.data.frame(as.table(corr))
  #remove the NA values from above 
  corr <- na.omit(corr) 
  #select significant values  
  corr <- subset(corr, abs(Freq) > sig) 
  #sort by highest correlation
  corr <- corr[order(-abs(corr$Freq)),] 
  #print table
  print(corr)
  #turn corr back into matrix in order to plot with corrplot
  mtx_corr <- reshape2::acast(corr, Var1 ~ Var2, value.var = "Freq")
  
  #plot correlations visually
  corrplot(mtx_corr, is.corr = FALSE, tl.col = "black", na.label = " ")
}
corr_simple()


# remove summary variables and variables that are highly correlated 

mlm_bl_analysis = subset(mlm_bl_analysis, select = -c(cSOPSPositive, cSOPSNegative, cSOPSGeneral, cSOPSDisorganization, cSOPSTotal, 
                                                      C_CDSTOTAL, GFS_current, GFS_lowest_pastyr, GFR_current, GFR_lowest_pastyr, GFS_highest_pastyr, GFR_highest_pastyr))
mlm_long_analysis = subset(mlm_long_analysis, select = -c(cSOPSPositive, cSOPSNegative, cSOPSGeneral, cSOPSDisorganization, cSOPSTotal, 
                                                          C_CDSTOTAL, GFS_current, GFS_lowest_pastyr, GFR_current, GFR_lowest_pastyr, GFS_highest_pastyr, GFR_highest_pastyr))

mlm_bl_analysis$converters = factor(mlm_bl_analysis$converters)
mlm_long_analysis$converters = factor(mlm_long_analysis$converters)
mlm_bl_analysis$VisitLabel = factor(mlm_bl_analysis$VisitLabel)
mlm_long_analysis$VisitLabel = factor(mlm_long_analysis$VisitLabel)

## prep data to run MLM for feature selection 
clin_var_names = colnames(mlm_long_analysis)[c(19:81,101:106)] # only include alcohol, THC, cigarettes
no_variables = length(clin_var_names)

fitlist = list()
fitlist = as.list(1:no_variables)
names(fitlist) = clin_var_names
pvalues = list()


for (i in clin_var_names) {
  # for loop to perform MLM on each possible predictor variable 
  tmp = mlm_long_analysis[, c(i, "VisitLabel", "converters", "NAPLS_ID")]
  fml = as.formula(paste(i, "~", paste(c("VisitLabel", "converters"), collapse = "*")))
  
  
  fit = lme(fml, data = mlm_long_analysis, 
            random = ~1 | NAPLS_ID, na.action = na.omit)
  # print(anova(fit)[[4,4]])
  fitlist[[i]] = lme(fml, data = mlm_long_analysis, 
                         random = ~1 | NAPLS_ID, na.action = na.omit)
}

for (i in fitlist) {
  # for loop to extract p-values of interaction term ONLY 
  pvalues = rbind(pvalues, anova(i)[[4,4]])
}

fit_test = lme(D1_SOPS ~ VisitLabel*converters, data = mlm_long_analysis,
               random = ~1 | NAPLS_ID, na.action = na.omit)

summary(fit_test)

# combine variables and p-values 
anovalist = data.frame(unlist(clin_var_names), unlist(pvalues))
anovalist$unlist.clin_var_names.

### Look at which variables are significant from this filtering process

# filter to just significant at 0.05 level 
select_sig_05 = filter(anovalist, unlist.pvalues. < 0.05)
select_sig_05$unlist.clin_var_names.
View(select_sig_05)

### FDR correction
fdr_list = select_sig_05 %>% 
  arrange(unlist.pvalues.)

fdr_list$fdr_threshold = 0

length = nrow(fdr_list)

for (i in 1:length) {
  fdr_list$fdr_threshold[i] = 0.05*(i/length)
}

fdr_list$include = ifelse(fdr_list$unlist.pvalues. < fdr_list$fdr_threshold, 1, 0)

fdr_filt = filter(fdr_list, include == 1)
fdr_filt

## filter data to variables that survived FDR correction in feature selection process:
jm_long_analysis = subset(mlm_long_analysis, select = c(NAPLS_ID, converters, days_from_baseline_until_convers, VisitLabel,
                                                        GFS_decline, psqinh_11, D1_SOPS, P3_SOPS))

write.csv(jm_long_analysis, "/Users/maw/Documents/Documents/Yale/Cannon Lab/Joint Modeling Risk Pred/Conversion/Data/jm_long_analysis_11.10.22.csv")


