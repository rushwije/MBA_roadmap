##################################################################################
#  Methodology for systematic identification of multiple bias analysis in causal 
#  inference : simulation of data  
#  Function sim.gen: for simulating data mimicking HealthNuts 
#  Author: Rushani Wijesuriya 29th of May 2024                            #
##################################################################################


library('fastDummies')
library(dplyr)

# #set working directory 
setwd("C:/Users/rushani.wijesuriya/OneDrive - Murdoch Children's Research Institute/QBA-methods paper/Simulation study")

#set working directory to where parameters.csv file is saved

#------------------------------------------------------------------------------#
#                          Function to generate data                           #  
#------------------------------------------------------------------------------#


#------ key inputs/parameters


#1.n <- sample size
#2.nsim <- number of simulations
#3.seed <- random seed
#4.p_U <- prevalance of the unmeasured confounder
#5.OR_AU <- Key bias parameter 1 (conditional association between exposure
#           and unmeasured confounder- OR)
#6.OR_YU <- Key bias parameter 2: (conditional association between outcome
#            and unmeasured confounder- OR) 
#7.OR_YAE <- #Key bias parameter 3:interaction between selection and exposure/effect
            # measure modification of the association due to S (will need to tinker)
            #- OR
#8.p_E <- non-selection probability into the study sample 
#9.sens_x <- sensitivity of exposure
#10.spec_x <- specificity of exposure
#11. sens_Y <- sensitivity of outcome
#12. spec_Y <- specificity of outcome
#11.OR_AstarA <- defined using sens and spec
#12. OR_YstarY <- defined using sens and spec
#13. OR_YRy <- probability of outcome missingness 
#14. alpha <- exposure coefficient (OR) in the outcome generation model 


#-----Outputs

#$result_text <-number of datasets omitted due to perfect seperation
#$dataframes <- Simulated data post selection 
#$dataframes_S <- Simulated data pre selection 


sim.gen <- function(n,nsim,seed,p_U,OR_AU,OR_YU,OR_YAE,p_E,sens_X,spec_X,sens_Y,spec_Y,OR_YRy,alpha){

set.seed(seed)
  
#import parameter values
param <- read.csv("parameters.csv")
param <- param[-1]
  
#function to invert logit 
  expit <- function(p) {
    exp(p) / (1 + exp(p))
  }

data.list <- list()
data.s.list <- list()

i=1
sep <- 0 #set a counter to capture datasets with perfect separation


while (length(Filter(Negate(is.null), data.list) )<nsim){

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


D$ccare <- sample(c(0,1),n,replace=TRUE,prob=c(param[param$para_name == 'CC_no', 2],
                                               param[param$para_name == 'CC_yes', 2]))

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


D$smokingii <- rbinom(n, 1, expit(param[param$para_name == 'smok_b0', 2] + 
                                    param[param$para_name == 'smok_b1', 2]*D$SEIFA_quint_1+
                                    param[param$para_name == 'smok_b2', 2]*D$SEIFA_quint_2+
                                    param[param$para_name == 'smok_b3', 2]*D$SEIFA_quint_3+
                                    param[param$para_name == 'smok_b4', 2]*D$SEIFA_quint_4 ))

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
                                      param[param$para_name == 'LBW.smok_b8', 2]*D$smokingii+
                                      param[param$para_name == 'LBW.msmok_b9', 2]*D$smoke_mother_pregnancyii))


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
                                  param[param$para_name == 'exp.ccare_b11', 2]*D$ccare+
                                  param[param$para_name == 'exp.dmode_b12', 2]*D$mode_delivery+
                                  param[param$para_name == 'exp.gest_b13', 2]*D$pre_term+
                                  param[param$para_name == 'exp.mage_b14', 2]*D$mother_age_atbirth+
                                  param[param$para_name == 'exp.SES2_b15', 2]*D$SEIFA_quint_1+
                                  param[param$para_name == 'exp.SES3_b16', 2]*D$SEIFA_quint_2+
                                  param[param$para_name == 'exp.SES4_b17', 2]*D$SEIFA_quint_3+
                                  param[param$para_name == 'exp.SES5_b18', 2]*D$SEIFA_quint_4+
                                  param[param$para_name == 'exp.smok_b19', 2]*D$smokingii+
                                  param[param$para_name == 'exp.msmok_b20', 2]*D$smoke_mother_pregnancyii+
                                  log(OR_AU)*D$U_gh
))


#----------Step 8: Generate selection indicator for type II selection (S=1 selected to the sample/can speak english)

D$E <- rbinom(n, 1,p_E)
#sample(c(0,1),n,replace=TRUE,prob=c(0.09,0.91)) 

#----------Step 9: Generate true outcome: asthma at age 6 

exp_coeff <- log(alpha)   #to be changed to achieve 80% power, see R script "4. Checking power.R"
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
                                                param[param$para_name == 'out.ccare_b12', 2]*D$ccare+
                                                param[param$para_name == 'out.dmode_b13', 2]*D$mode_delivery+
                                                param[param$para_name == 'out.gest_b14', 2]*D$pre_term+
                                                param[param$para_name == 'out.mage_b15', 2]*D$mother_age_atbirth+
                                                param[param$para_name == 'out.SES2_b16', 2]*D$SEIFA_quint_1+
                                                param[param$para_name == 'out.SES3_b17', 2]*D$SEIFA_quint_2+
                                                param[param$para_name == 'out.SES4_b18', 2]*D$SEIFA_quint_3+
                                                param[param$para_name == 'out.SES5_b19', 2]*D$SEIFA_quint_4+
                                                param[param$para_name == 'out.smok_b20', 2]*D$smokingii+
                                                param[param$para_name == 'out.msmok_b21', 2]*D$smoke_mother_pregnancyii
                                              
))


#----------Step 10: Generate misclassified exposure
# sens <- 0.82
# spec <- 0.93
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
                                  param[param$para_name == 'exp.ccare_b11', 2]*D$ccare+
                                  param[param$para_name == 'exp.dmode_b12', 2]*D$mode_delivery+
                                  param[param$para_name == 'exp.gest_b13', 2]*D$pre_term+
                                  param[param$para_name == 'exp.mage_b14', 2]*D$mother_age_atbirth+
                                  param[param$para_name == 'exp.SES2_b15', 2]*D$SEIFA_quint_1+
                                  param[param$para_name == 'exp.SES3_b16', 2]*D$SEIFA_quint_2+
                                  param[param$para_name == 'exp.SES4_b17', 2]*D$SEIFA_quint_3+
                                  param[param$para_name == 'exp.SES5_b18', 2]*D$SEIFA_quint_4+
                                  param[param$para_name == 'exp.smok_b19', 2]*D$smokingii+
                                  param[param$para_name == 'exp.msmok_b20', 2]*D$smoke_mother_pregnancyii+
                                  log(OR_AstarA)*D$bf_true
))

#----------Step 11: Generate miss classified outcome

# sens_Y <- 0.80
# spec_Y <- 0.97
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
                                    param[param$para_name == 'out.ccare_b12', 2]*D$ccare+
                                    param[param$para_name == 'out.dmode_b13', 2]*D$mode_delivery+
                                    param[param$para_name == 'out.gest_b14', 2]*D$pre_term+
                                    param[param$para_name == 'out.mage_b15', 2]*D$mother_age_atbirth+
                                    param[param$para_name == 'out.SES2_b16', 2]*D$SEIFA_quint_1+
                                    param[param$para_name == 'out.SES3_b17', 2]*D$SEIFA_quint_2+
                                    param[param$para_name == 'out.SES4_b18', 2]*D$SEIFA_quint_3+
                                    param[param$para_name == 'out.SES5_b19', 2]*D$SEIFA_quint_4+
                                    param[param$para_name == 'out.smok_b20', 2]*D$smokingii+
                                    param[param$para_name == 'out.msmok_b21', 2]*D$smoke_mother_pregnancyii+
                                    log(OR_YstarY)*D$asthma_dx_current_6y
                                  
))

#----------Step 12: Generate missingeness in outcome

param[param$para_name == 'Ry_b0', 2] <- -2.0
# phi_6 <- log(1.20)
D$Ry <- rbinom(n, 1, expit(   param[param$para_name == 'Ry_b0', 2] + 
                                param[param$para_name == 'Ry.exp_b1', 2]* D$bf_true+
                                param[param$para_name == 'Ry.sex_b2', 2]*D$sexii+
                                param[param$para_name == 'Ry.sibs1_b3', 2]*D$number_sibs_1+
                                param[param$para_name == 'Ry.sibs2_b4', 2]*D$number_sibs_2+
                                param[param$para_name == 'Ry.sibs3_b5', 2]*D$number_sibs_3+
                                param[param$para_name == 'Ry.LBW_b6', 2]*D$LBW_infants+
                                param[param$para_name == 'Ry.cobasi_b7', 2]*D$parent_cob_3cat_1+
                                param[param$para_name == 'Ry.cobotr_b8', 2]*D$parent_cob_3cat_2+
                                param[param$para_name == 'Ry.msath_b9', 2]*D$mother_asthmaii+
                                param[param$para_name == 'Ry.fsath_b10', 2]*D$father_asthmaii+
                                param[param$para_name == 'Ry.fhalle_b11', 2]*D$family_history_allergyii+
                                param[param$para_name == 'Ry.ccare_b12', 2]*D$ccare+
                                param[param$para_name == 'Ry.dmode_b13', 2]*D$mode_delivery+
                                param[param$para_name == 'Ry.gest_b14', 2]*D$pre_term+
                                param[param$para_name == 'Ry.mage_b15', 2]*D$mother_age_atbirth+
                                param[param$para_name == 'Ry.SES2_b16', 2]*D$SEIFA_quint_1+
                                param[param$para_name == 'Ry.SES3_b17', 2]*D$SEIFA_quint_2+
                                param[param$para_name == 'Ry.SES4_b18', 2]*D$SEIFA_quint_3+
                                param[param$para_name == 'Ry.SES5_b19', 2]*D$SEIFA_quint_4+
                                param[param$para_name == 'Ry.smok_b20', 2]*D$smokingii+
                                param[param$para_name == 'Ry.msmok_b21', 2]*D$smoke_mother_pregnancyii+
                                log(OR_YRy)*D$asthma_dx_current_6y
                              
))



#set Y (mismeasured) to be missing (save original) 
D$out_star_comp <-D$out_star 
D$out_star <- ifelse(D$Ry==1,D$out_star,NA)

#filter study sample
D_E <- D #first save original
D <- D %>% filter(E==1)

#set Y missing
# D$asthma_dx_current_6y <- ifelse(D$Ry==1,D$asthma_dx_current_6y,NA)

#remove unwanted variables (inter, S, Ry, true values)

D$inter <- NULL
D$error <- NULL

#check for separation issues 

p<- table(D$out_star,D$bf_star)
if(p[1,1]>0 & p[1,2]>0 & p[2,1]>0 & p[2,2]>0){

#save the data (in a list)
data.list[[i]] <- D
data.s.list[[i]] <- D_E

} else {sep <- sep+1}

i <- i+1

}

#return(data.list)
X<- paste("number of data sets omitted due to perfect seperation=",sep)

return(list(result_text = X, dataframes = Filter(Negate(is.null), data.list), dataframes_E=Filter(Negate(is.null), data.s.list)))



}

# check
# results <- sim.gen(n=2000,nsim=1000,seed=13082023,p_U=0.1,phi_1=log(0.6),phi_2=log(1.4),phi_3=log(0.7),p_S=0.8,
#         sens_X=0.82,spec_X=0.93, sens_Y=0.8,spec_Y=0.97,phi_6=log(1.20))
# 
# print(results)

