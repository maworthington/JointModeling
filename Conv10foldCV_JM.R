library(mice)
library(magrittr) 
library(JM)
library(pec)
library(rms)
library(caret)
library(tidyverse)
library(ggplot2)
library(modelr)
library(purrr)
library(parallel)
library(nlme)

## This script performs a 10-fold cross-validation on joint models developed 
## with the JM package. Feature selection was performed in a previous 
## analysis step and individual joint models are built for each 
## identified variable of interest. 


## Data importation and cleaning 

# import imputed data for joint modeling in long format
jm_imp_long_conv = read.csv("/Users/maw/Documents/Documents/Yale/Cannon Lab/Joint Modeling Risk Pred/Conversion/Data/JM_longdata_conversion_imputed.csv")

# convert variables to numeric and factor structure as appropriate 
jm_imp_long_conv[,c(4:14,16)] = sapply(jm_imp_long_conv[,c(4:14,16)], as.numeric)
jm_imp_long_conv$converters = as.factor(jm_imp_long_conv$converters)
jm_imp_long_conv$genetic = as.factor(jm_imp_long_conv$genetic)

# sort data by ID and observation time 
jm_imp_long_conv = jm_imp_long_conv %>% 
  arrange(NAPLS_ID, month)

# create wide version of the data set from the cleaned long version 
jm_imp_wide_conv = pivot_wider(jm_imp_long_conv, 
                          names_from = "month",
                          names_sep = ".",
                          values_from = c(5:8))

# create data to store performance results during CV process:
# values to be stored for each fold include AUC, mean risk score cutoff, sensitivity,
# specificity and balanced accuracy 
P3.results = data.frame(AUC = double(), 
                        mean_cutoff = double(), 
                        Sens = double(), 
                        Spec = double(), 
                        BAC = double(),
                        loglikelihood = double(),
                        AIC = double(),
                        BIC = double())

D1.results = data.frame(AUC = double(), 
                        mean_cutoff = double(), 
                        Sens = double(), 
                        Spec = double(), 
                        BAC = double(),
                        loglikelihood = double(),
                        AIC = double(),
                        BIC = double())

GFSl.results = data.frame(AUC = double(), 
                          mean_cutoff = double(), 
                          Sens = double(), 
                          Spec = double(), 
                          BAC = double(),
                          loglikelihood = double(),
                          AIC = double(),
                          BIC = double())

psqinh.results = data.frame(AUC = double(), 
                            mean_cutoff = double(), 
                            Sens = double(), 
                            Spec = double(), 
                            BAC = double(),
                            loglikelihood = double(),
                            AIC = double(),
                            BIC = double())

bl.results = data.frame(AUC = double(), 
                        mean_cutoff = double(), 
                        Sens = double(), 
                        Spec = double(), 
                        BAC = double())


## prep data for 10-fold cross validation by splitting into 10 folds 
set.seed(NULL)
set.seed(1234)

# create 10 folds
folds = createDataPartition(jm_imp_wide_conv$converters, p = 0.7, times = 10, list = FALSE)

## Perform 10 fold cross validation step 

for(i in 1:10) {
  # Segment into testing/training
  trainingDataWideConv = jm_imp_wide_conv[folds[,i], ]
  testingDataWideConv = jm_imp_wide_conv[-folds[,i], ]
  
  # match IDs for wide and long format training/testing data
  trainingDataConv = filter(jm_imp_long_conv, jm_imp_long_conv$NAPLS_ID %in% trainingDataWideConv$NAPLS_ID)
  testingDataConv = filter(jm_imp_long_conv, jm_imp_long_conv$NAPLS_ID %in% testingDataWideConv$NAPLS_ID)
  
  # Create BL data
  trainingDataConv.BL = filter(trainingDataConv, month == 0)
  testingDataConv.BL = filter(testingDataConv, month == 0)
  
  # arrange all data by ID and month 
  trainingDataConv = trainingDataConv %>%
    mutate(NAPLS_ID = as.numeric(NAPLS_ID)) %>%
    arrange(NAPLS_ID, month)
  trainingDataConv.BL = trainingDataConv.BL %>%
    mutate(NAPLS_ID = as.numeric(NAPLS_ID)) %>%
    arrange(NAPLS_ID, month)
  
  testingDataConv = testingDataConv %>%
    mutate(NAPLS_ID = as.numeric(NAPLS_ID)) %>%
    arrange(NAPLS_ID, month)
  testingDataConv.BL = testingDataConv.BL %>%
    mutate(NAPLS_ID = as.numeric(NAPLS_ID)) %>%
    arrange(NAPLS_ID, month)
  
  # fit LME models 
  P3.lme = lme(P3_SOPS ~ converters*month,
               random = ~ month | NAPLS_ID, data = trainingDataConv,
               control = lmeControl(opt = 'optim'))
  D1.lme = lme(D1_SOPS ~ converters*month,
               random = ~ month | NAPLS_ID, data = trainingDataConv,
               control = lmeControl(opt = 'optim'))
  GFSl.lme = lme(GFS_decline ~ converters*month,
                 random = ~ month | NAPLS_ID, data = trainingDataConv,
                 control = lmeControl(opt = 'optim'))
  psqinh.lme = lme(psqinh_11 ~ converters*month,
                 random = ~ month | NAPLS_ID, data = trainingDataConv,
                 control = lmeControl(opt = 'optim'))
  
  # fit baseline CPH model 
  full.cph = coxph(Surv(months_to_conversion, as.factor(converters) == 1) ~ demo_age_ym + P1P2 + hvlttotal + bacsraw + 
                     GFS_change_pastyr + as.factor(genetic) + clife_events_total_Undesirability + totaltraumas + 
                     cluster(NAPLS_ID),
                   data = trainingDataConv.BL, model = TRUE, x = TRUE)
  
  P3.cph = coxph(Surv(months_to_conversion, as.factor(converters) == 1) ~ demo_age_ym + P1P2 + hvlttotal + bacsraw + 
                     GFS_change_pastyr + as.factor(genetic) + clife_events_total_Undesirability + totaltraumas + 
                     P3_SOPS + cluster(NAPLS_ID),
                   data = trainingDataConv.BL, model = TRUE, x = TRUE)
  
  D1.cph = coxph(Surv(months_to_conversion, as.factor(converters) == 1) ~ demo_age_ym + P1P2 + hvlttotal + bacsraw + 
                   GFS_change_pastyr + as.factor(genetic) + clife_events_total_Undesirability + totaltraumas + 
                   D1_SOPS + cluster(NAPLS_ID),
                 data = trainingDataConv.BL, model = TRUE, x = TRUE)
  
  GFSl.cph = coxph(Surv(months_to_conversion, as.factor(converters) == 1) ~ demo_age_ym + P1P2 + hvlttotal + bacsraw + 
                   GFS_change_pastyr + as.factor(genetic) + clife_events_total_Undesirability + totaltraumas + 
                   cluster(NAPLS_ID),
                 data = trainingDataConv.BL, model = TRUE, x = TRUE)
  
  psqinh.cph = coxph(Surv(months_to_conversion, as.factor(converters) == 1) ~ demo_age_ym + P1P2 + hvlttotal + bacsraw + 
                   GFS_change_pastyr + as.factor(genetic) + clife_events_total_Undesirability + totaltraumas + 
                   psqinh_11 + cluster(NAPLS_ID),
                 data = trainingDataConv.BL, model = TRUE, x = TRUE)
  
  full.cph.2 = cph(Surv(months_to_conversion, as.factor(converters) == 1) ~ demo_age_ym + P1P2 + hvlttotal + bacsraw + 
                       GFS_change_pastyr + genetic + clife_events_total_Undesirability + totaltraumas,
                     data = trainingDataConv.BL, model = TRUE, x = TRUE, y = TRUE, surv = TRUE)
  
  # fit joint models
  P3.joint = jointModel(P3.lme, P3.cph, timeVar = "month")
  D1.joint = jointModel(D1.lme, D1.cph, timeVar = "month")
  GFSl.joint = jointModel(GFSl.lme, GFSl.cph, timeVar = "month")
  psqinh.joint = jointModel(psqinh.lme, psqinh.cph, timeVar = "month")
  
  # make survival predictions from joint models
  P3.survfit = survfitJM(P3.joint, newdata = testingDataConv, type = "SurvProb", 
                         simulate = T, idVar = "NAPLS_ID")
  D1.survfit = survfitJM(D1.joint, newdata = testingDataConv, type = "SurvProb", 
                         simulate = T, idVar = "NAPLS_ID")
  GFSl.survfit = survfitJM(GFSl.joint, newdata = testingDataConv, type = "SurvProb", 
                           simulate = T, idVar = "NAPLS_ID")
  psqinh.survfit = survfitJM(psqinh.joint, newdata = testingDataConv, type = "SurvProb", 
                             simulate = T, idVar = "NAPLS_ID")
  bl.surv = predictSurvProb(full.cph.2, newdata = testingDataConv.BL, times = 2*12)
  
  # test predictions from joint models
  # create lists to store intermediate performance results for each model
  risk.P3 = list()
  risk.D1 = list()
  risk.GFSl = list()
  risk.psqinh = list()
  risk.bl = list()
  
  preds = data.frame(NAPLS_ID = testingDataConv.BL$NAPLS_ID)
  
  # store predicted risk of conversion/survival for each model
  for (j in 1:nrow(testingDataConv.BL)) {
    risk.P3[j] = list(P3.survfit$summaries[[j]][24,2])
    risk.D1[j] = list(D1.survfit$summaries[[j]][24,2])
    risk.GFSl[j] = list(GFSl.survfit$summaries[[j]][24,2])
    risk.psqinh[j] = list(psqinh.survfit$summaries[[j]][24,2])
    risk.bl[j] = list(bl.surv)
  }
  
  # match individual predicted risk scores with ID 
  preds = preds %>%
    mutate(P3.risk_surv = as.numeric(risk.P3),
           P3.risk_conv = 1 - as.numeric(risk.P3),
           D1.risk_surv = as.numeric(risk.D1),
           D1.risk_conv = 1 - as.numeric(risk.D1),
           GFSl.risk_surv = as.numeric(risk.GFSl),
           GFSl.risk_conv = 1 - as.numeric(risk.GFSl),
           psqinh.risk_surv = as.numeric(risk.psqinh),
           psqinh.risk_conv = 1 - as.numeric(risk.psqinh),
           bl.risk_surv = as.numeric(bl.surv),
           bl.risk_conv = 1 - as.numeric(bl.surv))
  
  converters = testingDataConv.BL[,c(1,2)]
  predictions_all = merge(converters, preds, by = "NAPLS_ID")
  predictions_all$converters = factor(predictions_all$converters)
  
  # store performance metrics in prediction conversion from each fold for each model 
  P3.roc.jm = pROC::roc(predictions_all$converters ~ predictions_all$P3.risk_conv)
  P3.results[i, 1] = pROC::auc(P3.roc.jm)
  P3.results[i, 2] = mean(predictions_all$P3.risk_conv)
  P3.cutoff = ifelse(predictions_all$P3.risk_conv >= mean(predictions_all$P3.risk_conv), 1, 0)
  P3.cfmtx = confusionMatrix(as.factor(P3.cutoff), predictions_all$converters, positive = "1")
  P3.results[i, 3] = P3.cfmtx$byClass[1]
  P3.results[i, 4] = P3.cfmtx$byClass[2]
  P3.results[i, 5] = P3.cfmtx$byClass[11]
  P3.results[i, 6] = summary(P3.joint)[5]
  P3.results[i, 7] = summary(P3.joint)[6]
  P3.results[i, 8] = summary(P3.joint)[7]
 
  D1.roc.jm = pROC::roc(predictions_all$converters ~ predictions_all$D1.risk_conv)
  D1.results[i, 1] = pROC::auc(D1.roc.jm)
  D1.results[i, 2] = mean(predictions_all$D1.risk_conv)
  D1.cutoff = ifelse(predictions_all$D1.risk_conv >= mean(predictions_all$D1.risk_conv), 1, 0)
  D1.cfmtx = confusionMatrix(as.factor(D1.cutoff), predictions_all$converters, positive = "1")
  D1.results[i, 3] = D1.cfmtx$byClass[1]
  D1.results[i, 4] = D1.cfmtx$byClass[2]
  D1.results[i, 5] = D1.cfmtx$byClass[11]
  D1.results[i, 6] = summary(D1.joint)[5]
  D1.results[i, 7] = summary(D1.joint)[6]
  D1.results[i, 8] = summary(D1.joint)[7]
  
  GFSl.roc.jm = pROC::roc(predictions_all$converters ~ predictions_all$GFSl.risk_conv)
  GFSl.results[i, 1] = pROC::auc(GFSl.roc.jm)
  GFSl.results[i, 2] = mean(predictions_all$GFSl.risk_conv)
  GFSl.cutoff = ifelse(predictions_all$GFSl.risk_conv >= mean(predictions_all$GFSl.risk_conv), 1, 0)
  GFSl.cfmtx = confusionMatrix(as.factor(GFSl.cutoff), predictions_all$converters, positive = "1")
  GFSl.results[i, 3] = GFSl.cfmtx$byClass[1]
  GFSl.results[i, 4] = GFSl.cfmtx$byClass[2]
  GFSl.results[i, 5] = GFSl.cfmtx$byClass[11]
  GFSl.results[i, 6] = summary(GFSl.joint)[5]
  GFSl.results[i, 7] = summary(GFSl.joint)[6]
  GFSl.results[i, 8] = summary(GFSl.joint)[7]
  
  psqinh.roc.jm = pROC::roc(predictions_all$converters ~ predictions_all$psqinh.risk_conv)
  psqinh.results[i, 1] = pROC::auc(psqinh.roc.jm)
  psqinh.results[i, 2] = mean(predictions_all$psqinh.risk_conv)
  psqinh.cutoff = ifelse(predictions_all$psqinh.risk_conv >= mean(predictions_all$psqinh.risk_conv), 1, 0)
  psqinh.cfmtx = confusionMatrix(as.factor(psqinh.cutoff), predictions_all$converters, positive = "1")
  psqinh.results[i, 3] = psqinh.cfmtx$byClass[1]
  psqinh.results[i, 4] = psqinh.cfmtx$byClass[2]
  psqinh.results[i, 5] = psqinh.cfmtx$byClass[11]
  psqinh.results[i, 6] = summary(psqinh.joint)[5]
  psqinh.results[i, 7] = summary(psqinh.joint)[6]
  psqinh.results[i, 8] = summary(psqinh.joint)[7]
 
  bl.roc = pROC::roc(predictions_all$converters ~ predictions_all$bl.risk_conv)
  bl.results[i, 1] = pROC::auc(bl.roc)
  bl.results[i, 2] = mean(predictions_all$bl.risk_conv)
  bl.cutoff = ifelse(predictions_all$bl.risk_conv >= mean(predictions_all$bl.risk_conv), 1, 0)
  bl.cfmtx = confusionMatrix(as.factor(bl.cutoff), predictions_all$converters, positive = "1")
  bl.results[i, 3] = bl.cfmtx$byClass[1]
  bl.results[i, 4] = bl.cfmtx$byClass[2]
  bl.results[i, 5] = bl.cfmtx$byClass[11]
}
