##################################################################################
#  Methodology for systematic identification of multiple bias analysis in causal 
#  inference : Checking power and setting effect sizes - for the realistic case
#  Author: Rushani Wijesuriya 15th of August 2023                          #
##################################################################################

rm(list=ls())


library('fastDummies')
source("gcomp_boot_functions.R")
library(boot)

#import parameter values
param <- read.csv("parameters.csv")
param <- param[-1]

expit <- function(p) {
  exp(p) / (1 + exp(p))
}


#realistic
n <- 2000
p_U <- 0.10
OR_AU <- 0.60
OR_YU <-1.30 
OR_YAE <- 0.7 
p_E <- 0.85
sens_X <- 0.90
spec_X <- 0.84
sens_Y <- 0.83
spec_Y <- 0.90
OR_YRy <- 1.2

#enhanced 
n <- 2000
p_U <- 0.20
OR_AU <- 0.20
OR_YU <-1.60 
OR_YAE <- 0.4 
p_E <- 0.70
sens_X <- 0.80
spec_X <- 0.68
sens_Y <- 0.66
spec_Y <- 0.80
OR_YRy <- 1.4


p_val <- c()
p_val_RD <- c()

RD <- c()


set.seed(13082023)


for(i in 1:1000){

  
  print(i)

  D <- data.frame(id=c(1:n))
  
  #-------------------Observed confounders (Steps 1-5)
  
  
  #--------Step 1: Generate family history of asthma, family history of allergic disease,
  #parental ethnicity, number of older siblings, childcare attendance, child’s sex,
  #and SES. 
  D$family_history_allergyii <- sample(c(0,1),n,replace=TRUE, prob=c(param[param$para_name == 'Fx_al_no', 2],
                                                                     param[param$para_name == 'Fx_al_yes', 2]) )
  
  D$mother_asthmaii <- sample(c(0,1),n,replace=TRUE,prob=c(param[param$para_name == 'M_A_no', 2],
                                                           param[param$para_name == 'M_A_yes', 2]))
  
  
  D$father_asthmaii <- sample(c(0,1),n,replace=TRUE,prob=c(param[param$para_name == 'F_A_no', 2],
                                                           param[param$para_name == 'F_A_yes', 2]))
  
  D$parent_cob_3cat <- sample(c(0,1,2),n,replace=TRUE,prob=c(param[param$para_name == 'CoB_aus', 2],
                                                             param[param$para_name == 'CoB_asia', 2],
                                                             param[param$para_name == 'CoB_other', 2]))
  
  
  D$number_sibs <- sample(c(0,1,2,3),n,replace=TRUE,prob=c(param[param$para_name == 'sib_no', 2],
                                                           param[param$para_name == 'sib_1', 2],
                                                           param[param$para_name == 'sib_2', 2],
                                                           param[param$para_name == 'sib_3', 2]))
  
  
  D$SEIFA_quint <- sample(c(0,1,2,3,4),n,replace=TRUE,prob=c(param[param$para_name == 'SES_1', 2],
                                                             param[param$para_name == 'SES_2', 2],
                                                             param[param$para_name == 'SES_3', 2],
                                                             param[param$para_name == 'SES_4', 2],
                                                             param[param$para_name == 'SES_5', 2]))
  
  
  D$sexii <-sample(c(0,1),n,replace=TRUE,prob=c(param[param$para_name == 'sex_m', 2],
                                                param[param$para_name == 'sex_f', 2])) 
  
  
  #--------Step 2: Generate maternal age at childbirth, env smoking and maternal smoking
  
  #create dummy indicators for SES
  
  D=fastDummies::dummy_cols(D, select_columns =c("parent_cob_3cat","number_sibs","SEIFA_quint"))
  
  D$error <- rnorm(n,0,4.8)
  D$mother_age_atbirth <-param[param$para_name == 'mage_b0', 2] +  param[param$para_name == 'mage_b1', 2]*D$SEIFA_quint_1+
    param[param$para_name == 'mage_b2', 2]*D$SEIFA_quint_2+
    param[param$para_name == 'mage_b3', 2]*D$SEIFA_quint_3+
    param[param$para_name == 'mage_b4', 2]*D$SEIFA_quint_4 + D$error
  
  D$smoke_mother_pregnancyii <- rbinom(n, 1, expit(param[param$para_name == 'msmok_b0', 2] + 
                                                     param[param$para_name == 'msmok_b1', 2]*D$SEIFA_quint_1+
                                                     param[param$para_name == 'msmok_b2', 2]*D$SEIFA_quint_2+
                                                     param[param$para_name == 'msmok_b3', 2]*D$SEIFA_quint_3+
                                                     param[param$para_name == 'msmok_b4', 2]*D$SEIFA_quint_4 ))
  
  
  #--------Step 3: Generate gestational age
  D$pre_term <- rbinom(n, 1, expit(param[param$para_name == 'gest_b0', 2] + 
                                     param[param$para_name == 'gest.mage_b1', 2]*D$mother_age_atbirth+
                                     param[param$para_name == 'gest.SES2_b2', 2]*D$SEIFA_quint_1+
                                     param[param$para_name == 'gest.SES3_b3', 2]*D$SEIFA_quint_2+
                                     param[param$para_name == 'gest.SES4_b4', 2]*D$SEIFA_quint_3+
                                     param[param$para_name == 'gest.SES5_b5', 2]*D$SEIFA_quint_4+
                                     param[param$para_name == 'gest.smok_b6', 2]*D$smoke_mother_pregnancyii))
  
  #--------Step 4: Generate mode of delivery
  
  D$mode_delivery <- rbinom(n, 1, expit(param[param$para_name == 'dmode_b0', 2] + 
                                          param[param$para_name == 'dmode.mage_b1', 2]*D$mother_age_atbirth+
                                          param[param$para_name == 'dmode.SES2_b2', 2]*D$SEIFA_quint_1+
                                          param[param$para_name == 'dmode.SES3_b3', 2]*D$SEIFA_quint_2+
                                          param[param$para_name == 'dmode.SES4_b4', 2]*D$SEIFA_quint_3+
                                          param[param$para_name == 'dmode.SES5_b5', 2]*D$SEIFA_quint_4))
  
  
  
  #---------Step 5: Generate birth weight 
  
  D$LBW_infants <- rbinom(n, 1, expit(param[param$para_name == 'LBW_b0', 2] + 
                                        param[param$para_name == 'LBW.dmode_b1', 2]*D$mode_delivery+
                                        param[param$para_name == 'LBW.gest_b2', 2]*D$pre_term+
                                        param[param$para_name == 'LBW.mage_b3', 2]*D$mother_age_atbirth+
                                        param[param$para_name == 'LBW.SES2_b4', 2]*D$SEIFA_quint_1+
                                        param[param$para_name == 'LBW.SES3_b5', 2]*D$SEIFA_quint_2+
                                        param[param$para_name == 'LBW.SES4_b6', 2]*D$SEIFA_quint_3+
                                        param[param$para_name == 'LBW.SES5_b7', 2]*D$SEIFA_quint_4+
                                        param[param$para_name == 'LBW.msmok_b8', 2]*D$smoke_mother_pregnancyii))
  
  
  #----------Step 6: Generate unmeasured confounder: gestational hypertension
  
  #D$U_gh <- sample(c(0,1),n,replace=TRUE,prob=c(0.90,0.10)) #Check 
  
  D$U_gh <- rbinom(n, 1,p_U)
  
  
  #-----------Step 7: Generate exposure: breast feeding (true)
  #phi_1 <-      #Key bias parameter 1 (conditional association between exposure
  # and unknown confounder)
  D$bf_true <- rbinom(n, 1, expit(param[param$para_name == 'exp_b0', 2] + 
                                    param[param$para_name == 'exp.sex_b1', 2]*D$sexii+
                                    param[param$para_name == 'exp.sibs1_b2', 2]*D$number_sibs_1+
                                    param[param$para_name == 'exp.sibs2_b3', 2]*D$number_sibs_2+
                                    param[param$para_name == 'exp.sibs3_b4', 2]*D$number_sibs_3+
                                    param[param$para_name == 'exp.LBW_b5', 2]*D$LBW_infants+
                                    param[param$para_name == 'exp.cobasi_b6', 2]*D$parent_cob_3cat_1+
                                    param[param$para_name == 'exp.cobotr_b7', 2]*D$parent_cob_3cat_2+
                                    param[param$para_name == 'exp.msath_b8', 2]*D$mother_asthmaii+
                                    param[param$para_name == 'exp.fsath_b9', 2]*D$father_asthmaii+
                                    param[param$para_name == 'exp.fhalle_b10', 2]*D$family_history_allergyii+
                                    param[param$para_name == 'exp.dmode_b11', 2]*D$mode_delivery+
                                    param[param$para_name == 'exp.gest_b12', 2]*D$pre_term+
                                    param[param$para_name == 'exp.mage_b13', 2]*D$mother_age_atbirth+
                                    param[param$para_name == 'exp.SES2_b14', 2]*D$SEIFA_quint_1+
                                    param[param$para_name == 'exp.SES3_b15', 2]*D$SEIFA_quint_2+
                                    param[param$para_name == 'exp.SES4_b16', 2]*D$SEIFA_quint_3+
                                    param[param$para_name == 'exp.SES5_b17', 2]*D$SEIFA_quint_4+
                                    param[param$para_name == 'exp.msmok_b18', 2]*D$smoke_mother_pregnancyii+
                                    log(OR_AU)*D$U_gh
  ))
  
  
  #----------Step 8: Generate selection indicator for type II selection (S=1 selected to the sample/can speak english)
  
  D$E <- rbinom(n, 1,p_E)
  #sample(c(0,1),n,replace=TRUE,prob=c(0.09,0.91)) 
  
  #----------Step 9: Generate true outcome: asthma at age 6 
  
  exp_coeff <- log(0.95)   #to be changed to achieve 80% power (0.78 in data)
  E_coeff <-   log(1.2)     #association between S and asthma 

  #Generate interaction term 
  
  D$inter <- D$E*D$bf_true
  
  
  D$asthma_dx_current_6y <- rbinom(n, 1, expit( exp_coeff*D$bf_true+   #exposure
                                                  E_coeff*D$E+         #selection variable
                                                  log(OR_YU)*D$U_gh+       #unmeasured confounder
                                                  log(OR_YAE)*D$inter+   #interaction term
                                                  param[param$para_name == 'out_b0', 2] + 
                                                  param[param$para_name == 'out.sex_b2', 2]*D$sexii+
                                                  param[param$para_name == 'out.sibs1_b3', 2]*D$number_sibs_1+
                                                  param[param$para_name == 'out.sibs2_b4', 2]*D$number_sibs_2+
                                                  param[param$para_name == 'out.sibs3_b5', 2]*D$number_sibs_3+
                                                  param[param$para_name == 'out.LBW_b6', 2]*D$LBW_infants+
                                                  param[param$para_name == 'out.cobasi_b7', 2]*D$parent_cob_3cat_1+
                                                  param[param$para_name == 'out.cobotr_b8', 2]*D$parent_cob_3cat_2+
                                                  param[param$para_name == 'out.msath_b9', 2]*D$mother_asthmaii+
                                                  param[param$para_name == 'out.fsath_b10', 2]*D$father_asthmaii+
                                                  param[param$para_name == 'out.fhalle_b11', 2]*D$family_history_allergyii+
                                                  param[param$para_name == 'out.dmode_b12', 2]*D$mode_delivery+
                                                  param[param$para_name == 'out.gest_b13', 2]*D$pre_term+
                                                  param[param$para_name == 'out.mage_b14', 2]*D$mother_age_atbirth+
                                                  param[param$para_name == 'out.SES2_b15', 2]*D$SEIFA_quint_1+
                                                  param[param$para_name == 'out.SES3_b16', 2]*D$SEIFA_quint_2+
                                                  param[param$para_name == 'out.SES4_b17', 2]*D$SEIFA_quint_3+
                                                  param[param$para_name == 'out.SES5_b18', 2]*D$SEIFA_quint_4+
                                                  param[param$para_name == 'out.msmok_b19', 2]*D$smoke_mother_pregnancyii
                                                
  ))
  
  
  #----------Step 10: Generate misclassified exposure
  # sens <- 0.82
  # spec <- 0.93
  # table(D$bf_true)
  # 
  # TP <- table(D$bf_true)[2]*sens_X
  # FN <-  table(D$bf_true)[2]-TP
  # 
  # FP <- table(D$bf_true)[1]*spec_X
  # TN <- table(D$bf_true)[1]-FP
  # 
  p_A <- 0.94
  
  TP <- p_A*n*sens_X
  FN <-  (p_A*n)-TP
  
  FP <- (1-p_A)*n*spec_X
  TN <- ((1-p_A)*n)-FP
  
  OR_AstarA <-  (TP/FN)/(FP/TN)
  D$bf_star <- rbinom(n, 1, expit(param[param$para_name == 'exp_b0', 2] + 
                                    param[param$para_name == 'exp.sex_b1', 2]*D$sexii+
                                    param[param$para_name == 'exp.sibs1_b2', 2]*D$number_sibs_1+
                                    param[param$para_name == 'exp.sibs2_b3', 2]*D$number_sibs_2+
                                    param[param$para_name == 'exp.sibs3_b4', 2]*D$number_sibs_3+
                                    param[param$para_name == 'exp.LBW_b5', 2]*D$LBW_infants+
                                    param[param$para_name == 'exp.cobasi_b6', 2]*D$parent_cob_3cat_1+
                                    param[param$para_name == 'exp.cobotr_b7', 2]*D$parent_cob_3cat_2+
                                    param[param$para_name == 'exp.msath_b8', 2]*D$mother_asthmaii+
                                    param[param$para_name == 'exp.fsath_b9', 2]*D$father_asthmaii+
                                    param[param$para_name == 'exp.fhalle_b10', 2]*D$family_history_allergyii+
                                    param[param$para_name == 'exp.dmode_b11', 2]*D$mode_delivery+
                                    param[param$para_name == 'exp.gest_b12', 2]*D$pre_term+
                                    param[param$para_name == 'exp.mage_b13', 2]*D$mother_age_atbirth+
                                    param[param$para_name == 'exp.SES2_b14', 2]*D$SEIFA_quint_1+
                                    param[param$para_name == 'exp.SES3_b15', 2]*D$SEIFA_quint_2+
                                    param[param$para_name == 'exp.SES4_b16', 2]*D$SEIFA_quint_3+
                                    param[param$para_name == 'exp.SES5_b17', 2]*D$SEIFA_quint_4+
                                    param[param$para_name == 'exp.msmok_b18', 2]*D$smoke_mother_pregnancyii+
                                    log(OR_AstarA)*D$bf_true
  ))
  
  #----------Step 11: Generate miss classified outcome
  
  # sens_Y <- 0.80
  # spec_Y <- 0.97
  # table(D$asthma_dx_current_6y)
  # 
  # TP <- table(D$asthma_dx_current_6y)[2]*sens_Y
  # FN <-  table(D$asthma_dx_current_6y)[2]-TP
  # 
  # FP <- table(D$asthma_dx_current_6y)[1]*spec_Y
  # TN <- table(D$asthma_dx_current_6y)[1]-FP
  # 
  p_Y <- 0.10
  
  TP <- p_Y*n*sens_Y
  FN <-  (p_Y*n)-TP
  
  FP <- (1-p_Y)*n*spec_Y
  TN <- ((1-p_Y)*n)-FP
  
  OR_YstarY <-  (TP/FN)/(FP/TN)
  D$out_star <- rbinom(n, 1, expit( param[param$para_name == 'out_b0', 2] + 
                                      param[param$para_name == 'out.sex_b2', 2]*D$sexii+
                                      param[param$para_name == 'out.sibs1_b3', 2]*D$number_sibs_1+
                                      param[param$para_name == 'out.sibs2_b4', 2]*D$number_sibs_2+
                                      param[param$para_name == 'out.sibs3_b5', 2]*D$number_sibs_3+
                                      param[param$para_name == 'out.LBW_b6', 2]*D$LBW_infants+
                                      param[param$para_name == 'out.cobasi_b7', 2]*D$parent_cob_3cat_1+
                                      param[param$para_name == 'out.cobotr_b8', 2]*D$parent_cob_3cat_2+
                                      param[param$para_name == 'out.msath_b9', 2]*D$mother_asthmaii+
                                      param[param$para_name == 'out.fsath_b10', 2]*D$father_asthmaii+
                                      param[param$para_name == 'out.fhalle_b11', 2]*D$family_history_allergyii+
                                      param[param$para_name == 'out.dmode_b12', 2]*D$mode_delivery+
                                      param[param$para_name == 'out.gest_b13', 2]*D$pre_term+
                                      param[param$para_name == 'out.mage_b14', 2]*D$mother_age_atbirth+
                                      param[param$para_name == 'out.SES2_b15', 2]*D$SEIFA_quint_1+
                                      param[param$para_name == 'out.SES3_b16', 2]*D$SEIFA_quint_2+
                                      param[param$para_name == 'out.SES4_b17', 2]*D$SEIFA_quint_3+
                                      param[param$para_name == 'out.SES5_b18', 2]*D$SEIFA_quint_4+
                                      param[param$para_name == 'out.msmok_b19', 2]*D$smoke_mother_pregnancyii+
                                      log(OR_YstarY)*D$asthma_dx_current_6y
                                    
  ))
  
  
  model <-  glm(asthma_dx_current_6y~bf_true+U_gh+sexii+as.factor(number_sibs)+mother_age_atbirth+
                  LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                  family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                  smoke_mother_pregnancyii, data=D, family=binomial(logit))
  
 p_val[i] <- summary(model)[["coefficients"]][2,4]
 
 
 
 model_RD <- glm(asthma_dx_current_6y~bf_true+sexii+as.factor(number_sibs)+mother_age_atbirth+
                   LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                   family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                   smoke_mother_pregnancyii+U_gh,family = gaussian("identity"), data=D)
 
 p_val_RD[i] <- summary(model_RD)[["coefficients"]][2,4]
 RD[i] <- summary(model_RD)[["coefficients"]][2,1]
 
 
}


#check power
sum(p_val<0.05)/1000
sum(p_val_RD<0.05)/1000
mean(RD)

