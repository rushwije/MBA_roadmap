##################################################################################
#  Application of QBA techniques to the HealthNuts case study
#  Data used: HealthNuts data 
#  Author: Rushani Wijesuriya 
#  19th of November 2024                             
##################################################################################

library(tidyverse)
library(doParallel)
library(readstata13)
library(gtools)
library('fastDummies')
library(survey)
library(sandwich)


#function for inverse logit
expit <- function(p) {
  exp(p) / (1 + exp(p))
}



#function for drawing bias paramaeters from N(mean,sd)
normBP <- function(i,model){
  
  mu <- summary(model)$coef[i,1]
  sd <- summary(model)$coef[i,2]
  return(c(mu,sd))
}


#set working directory
#setwd("~/case study application")


#import data 
dat <- read.dta13("HN_BF_asthma_clean.dta", nonint.factors = TRUE,generate.factors = T)
dat <- dat %>% select("id" ,"asthma_dx_current_6y" ,"sexii" ,"mode_delivery",
                      "pre_term","birth_weight","any_brfeed_an","mother_age_atbirth",
                      "number_sibs","smoke_mother_pregnancyii",
                      "SEIFA_quint" ,"parent_cob_3cat","mother_asthmaii",
                      "father_asthmaii","family_history_allergyii",
                     "quest_type_h3","brfeed")


#-----------------preliminary data recoding


#generate birth weight variable 
dat$LBW_infants <- ifelse(dat$birth_weight>=2500,0,1)

#change all factor variables to correct type
cat.vars <- c("asthma_dx_current_6y","smoke_mother_pregnancyii",
              "SEIFA_quint","mother_asthmaii","father_asthmaii",
              "family_history_allergyii","LBW_infants" )

dat[cat.vars]=lapply(dat[cat.vars],factor)

dat$any_brfeed_an <- relevel(dat$any_brfeed_an,ref="no")

confounders <- c("sexii" ,"mode_delivery","pre_term","LBW_infants",
                 "mother_age_atbirth","smoke_mother_pregnancyii","parent_cob_3cat",
                 "SEIFA_quint" ,"mother_asthmaii","father_asthmaii","family_history_allergyii",
                 "number_sibs")

#prep data and generate indicator variables for the bias models 
dat <- dummy_cols(dat, select_columns = c('SEIFA_quint',"parent_cob_3cat","number_sibs"),
                  remove_first_dummy = TRUE)
dat <- rename(dat, parent_cob_3cat_1_or_both_asian = "parent_cob_3cat_1 or both asian",
              number_sibs_1_sibling="number_sibs_1 sibling", 
              number_sibs_2_siblings ="number_sibs_2 siblings",
              number_sibs_3_or_more_siblings="number_sibs_3 or more siblings")

cat.vars2 <- c("any_brfeed_an","sexii","mode_delivery","pre_term")

dat[cat.vars2]=lapply(dat[cat.vars2],as.numeric)
dat[cat.vars2] <- dat[cat.vars2]-1

#Set up Table 5 
Table5 <- data.frame(matrix(ncol = 5, nrow = 7))
names(Table5) <- c("Biases adjusted for","Risk difference","95%CI","Risk ratio","95%CI")

#-------------------------------------------------------------------------------
#  Primary analysis
#-------------------------------------------------------------------------------

dat_prim <- dat
outcome <- c("asthma_dx_current_6y")
exp_confounders <- c("any_brfeed_an","sexii" ,"mode_delivery","pre_term","LBW_infants",
                     "mother_age_atbirth","smoke_mother_pregnancyii","parent_cob_3cat",
                     "SEIFA_quint" ,"mother_asthmaii","father_asthmaii","family_history_allergyii",
                     "number_sibs") 

confounders <- c("sexii" ,"mode_delivery","pre_term","LBW_infants",
                 "mother_age_atbirth","smoke_mother_pregnancyii","parent_cob_3cat",
                 "SEIFA_quint" ,"mother_asthmaii","father_asthmaii","family_history_allergyii",
                 "number_sibs")


#Risk difference 
Fit2 <- glm(as.formula(paste(outcome, paste(exp_confounders, collapse=" + "), sep=" ~ ")),
            family = binomial("identity"), dat_prim)

summary(Fit2)$coef[2, 1] 
RD_CI <- round(c((summary(Fit2)$coef[2, 1] + summary(Fit2)$coef[2, 2] * qnorm(.025)), 
           summary(Fit2)$coef[2, 1] + summary(Fit2)$coef[2, 2] * qnorm(.975)),2)


#Risk ratio
dat_prim$asthma_dx_current_6y <- as.numeric(dat_prim$asthma_dx_current_6y)
Fit3 <- glm(as.formula(paste(outcome, paste(exp_confounders, collapse=" + "), sep=" ~ ")),
            family = poisson("log"), dat_prim)

exp(summary(Fit3)$coef[2, 1] )
RR_CI <- round(exp(c((summary(Fit3)$coef[2, 1] + summary(Fit3)$coef[2, 2] * qnorm(.025)), 
               summary(Fit3)$coef[2, 1] + summary(Fit3)$coef[2, 2] * qnorm(.975))),2)


Table5[1,] <- c("Primary analysis",round(summary(Fit2)$coef[2, 1],2),
                paste(RD_CI[1],RD_CI[2]),round(exp(summary(Fit3)$coef[2, 1]),2)
,paste(RR_CI[1],RR_CI[2]))

#-------------------------------------------------------------------------------
#1. Selection bias-Type I (Collider stratification bias)
#-------------------------------------------------------------------------------


set.seed(38034) #for replicability 

#----- Step1 : obtain  bias parameter estimates 

#Note here we estimate the bias parameters for the exposure and confounders 
#based on observed data and for the outcome based on expert opinion   

#1. Generate missingness indicator (missing=1,observed=0) and fit a logistic
#regression model to obtain bias parameters for bias models 
#NOTE: Alternativey can also obtain the marginal estimates for each variables and 
# dilute their magnitude 

dat$miss <- ifelse(dat$quest_type_h3=="full q",1,0)

outcome <- c("miss")
all.vars <- c("any_brfeed_an","sexii" ,"mode_delivery","pre_term","LBW_infants",
                     "mother_age_atbirth","smoke_mother_pregnancyii","parent_cob_3cat",
              "SEIFA_quint" ,"mother_asthmaii","father_asthmaii","family_history_allergyii",
                     "number_sibs")  

propmodel <- glm(as.formula(paste(outcome, paste(all.vars, collapse=" + "), sep=" ~ ")),
                 family=binomial(logit), data=dat)


mu_y <- log(1.2)  
sd_y <-  (log(1.5)-log(1.1))/(2*1.96)

hist(exp(rnorm(10000,mu_y,sd_y)))


#vector to save results 
coef_RD <- c()
coef_RR <- c()



#no of bootstrap samples
I <- 1000

#step 2: fit the bias model 
for(i in 1:I){
  
  print(i)
  bootstrap_indices <- sample(1:nrow(dat), nrow(dat), replace=TRUE)

  s1_0     <- rnorm(1,normBP(1,propmodel)[1],normBP(1,propmodel)[2])
  s1_x     <- rnorm(1, normBP(2,propmodel)[1],normBP(2,propmodel)[2])
  
  s1_y     <- rnorm(1,mean=mu_y,sd=sd_y)
  
  s1_c1    <-rnorm(1, normBP(3,propmodel)[1],normBP(3,propmodel)[2])
  s1_c2    <-rnorm(1, normBP(4,propmodel)[1],normBP(4,propmodel)[2])
  s1_c3    <-rnorm(1, normBP(5,propmodel)[1],normBP(5,propmodel)[2])
  s1_c4    <-rnorm(1, normBP(6,propmodel)[1],normBP(6,propmodel)[2])
  s1_c5    <-rnorm(1, normBP(7,propmodel)[1],normBP(7,propmodel)[2])
  s1_c6    <-rnorm(1, normBP(8,propmodel)[1],normBP(8,propmodel)[2])
  s1_c71    <-rnorm(1, normBP(9,propmodel)[1],normBP(9,propmodel)[2])
  s1_c72    <-rnorm(1, normBP(10,propmodel)[1],normBP(10,propmodel)[2])
  
  s1_c81    <-rnorm(1, normBP(11,propmodel)[1],normBP(11,propmodel)[2])
  s1_c82    <-rnorm(1, normBP(12,propmodel)[1],normBP(12,propmodel)[2])
  s1_c83    <-rnorm(1, normBP(13,propmodel)[1],normBP(13,propmodel)[2])
  s1_c84    <-rnorm(1, normBP(14,propmodel)[1],normBP(14,propmodel)[2])
  
  s1_c9    <-rnorm(1, normBP(15,propmodel)[1],normBP(15,propmodel)[2])
  s1_c10   <-rnorm(1, normBP(16,propmodel)[1],normBP(16,propmodel)[2])
  s1_c11   <-rnorm(1, normBP(17,propmodel)[1],normBP(17,propmodel)[2])
  
  s1_c121   <-rnorm(1, normBP(18,propmodel)[1],normBP(18,propmodel)[2])
  s1_c122   <-rnorm(1, normBP(19,propmodel)[1],normBP(19,propmodel)[2])
  s1_c123   <-rnorm(1, normBP(20,propmodel)[1],normBP(20,propmodel)[2])
  

  df <- dat[bootstrap_indices,]
  #code to numeric 
  df[cat.vars]=lapply(df[cat.vars],as.numeric)
  df[cat.vars]= df[cat.vars]-1
  

df <- df  %>% 
  mutate( pS = expit(s1_0 + s1_x * any_brfeed_an + s1_y * asthma_dx_current_6y + s1_c1 * sexii+s1_c2 * mode_delivery+
                       s1_c3 * pre_term+s1_c4 * LBW_infants+s1_c5 * mother_age_atbirth+ s1_c6 * smoke_mother_pregnancyii+
                       s1_c71 * parent_cob_3cat_1_or_both_asian + s1_c71*parent_cob_3cat_other+ s1_c81* SEIFA_quint_2+
                       s1_c82* SEIFA_quint_3+ s1_c83* SEIFA_quint_4+ s1_c84* SEIFA_quint_5 +s1_c9* mother_asthmaii+
                       s1_c10 * father_asthmaii+s1_c11 * family_history_allergyii+
                       s1_c121 * number_sibs_1_sibling +s1_c122 *number_sibs_2_siblings +s1_c122* number_sibs_3_or_more_siblings) )

cat.vars_an <- names(df)[c(3:7,9:20)]
df[cat.vars_an]=lapply(df[cat.vars_an],factor)

df <- df %>% filter(!is.na(pS))

outcome <- c("asthma_dx_current_6y")
exp_confounders <- c("any_brfeed_an","sexii","mode_delivery","pre_term","LBW_infants",
                     "mother_age_atbirth","smoke_mother_pregnancyii","parent_cob_3cat",
                     "SEIFA_quint" ,"mother_asthmaii","father_asthmaii","family_history_allergyii",
                     "number_sibs")

#-----------------risk difference
df$ipw <- 1/df$pS
df <- df %>% filter(!is.na(ipw))

final_RD <- svyglm(as.formula(paste(outcome, paste(exp_confounders, collapse=" + "), sep=" ~ ")), family=gaussian("identity"),
                   design=svydesign(id=~1, weights=~ipw, data=df))

coef_RD[i] <- coef(final_RD)[2]

#-------------------Risk ratio
final_RR <- svyglm(as.formula(paste(outcome, paste(exp_confounders, collapse=" + "), sep=" ~ ")), family=poisson("log"),
                   design=svydesign(id=~1, weights=~ipw, data=df))


coef_RR[i] <-coef(final_RR)[2]


}


RD <- round(median(coef_RD),2) #estimate
RD_CI <- round(quantile(coef_RD,c(.025, .975)),2) # 95% CI


RR <- round(exp(median(coef_RR)),2) # estimate
RR_CI <- round(exp(quantile(coef_RR,c(.025, .975))),2) # 95% CI 

Table5[3,] <- c("SB-collider stratification",RD,
                paste(RD_CI[1],RD_CI[2]),RR,paste(RR_CI[1],RR_CI[2]))



#----------------------------------------------------------
#2. Selection bias- Type II (generalizability bias)- consent 
#---------------------------------------------------------

set.seed(94782352) #for replicability 

#----- Step1 : obtain  bias parameter estimates 

S_0mu <- log(2.8)
S_0sd <- (log(4)-log(2.3))/(2*1.96)

S_xmu <- log(1.1)
S_xsd <- (log(1.3)-log(1.0))/(2*1.96)

S_ymu <- log(1.2)
S_ysd <- (log(1.3)-log(0.90))/(2*1.96)

# S_xymu <-log(1.3)
# S_xysd <- (log(1.6)-log(1.1))/(2*1.96)

S_xymu <-log(1.4)
S_xysd <- (log(1.7)-log(1.0))/(2*1.96)

hist(exp(rnorm(10000,S_xymu,S_xysd)))



#vector to save results 
coef_RD <- c()
coef_RR <- c()


#no of bootstrap samples
I <- 1000

#step 2: fit the bias model 

for(i in 1:I){
  
  print(i)
  bootstrap_indices <- sample(1:nrow(dat), nrow(dat), replace=TRUE)
  
  
  S_0 <- rnorm(1,S_0mu,S_0sd)
  S_x <- rnorm(1,S_xmu,S_xsd)
  S_y <- rnorm(1,S_ymu,S_ysd)
  S_xy <- rnorm(1,S_xymu,S_xysd)
  
  df <- dat[bootstrap_indices,]
  
  df[cat.vars]=lapply(df[cat.vars],as.numeric)
  df[cat.vars]= df[cat.vars]-1
  
  df <-  df %>% 
    mutate( pS = expit(S_0 +S_x  * any_brfeed_an + S_y * asthma_dx_current_6y + S_xy* asthma_dx_current_6y*any_brfeed_an )
    )
  
  cat.vars_an <- names(df)[c(3:7,9:17)]
  df[cat.vars_an]=lapply(df[cat.vars_an],factor)
  
  
  
  outcome <- c("asthma_dx_current_6y")
  exp_confounders <- c("any_brfeed_an","sexii","mode_delivery","pre_term","LBW_infants",
                       "mother_age_atbirth","smoke_mother_pregnancyii","parent_cob_3cat",
                       "SEIFA_quint" ,"mother_asthmaii","father_asthmaii","family_history_allergyii",
                       "number_sibs")   
  
  
  df$ipw <- 1/df$pS
  
  df <- df %>% filter(!is.na(ipw))
  final_RD <- svyglm(as.formula(paste(outcome, paste(exp_confounders, collapse=" + "), sep=" ~ ")), family=gaussian("identity"),
                     design=svydesign(id=~1, weights=~ipw, data=df))
  
  coef_RD[i] <- coef(final_RD)[2]
  
  #-------------------Risk ratio
  final_RR <- svyglm(as.formula(paste(outcome, paste(exp_confounders, collapse=" + "), sep=" ~ ")), family=poisson("log"),
                     design=svydesign(id=~1, weights=~ipw, data=df))
  
  coef_RR[i] <- coef(final_RR)[2]
  
}



RD <- round(median(coef_RD),2) #estimate
RD_CI <- round(quantile(coef_RD,c(.025, .975)),2) # 95% CI


RR <- round(exp(median(coef_RR)),2) # estimate
RR_CI <- round(exp(quantile(coef_RR,c(.025, .975))),2) # 95% CI 

Table5[4,] <- c("SB-generalizability-consent",RD,
                paste(RD_CI[1],RD_CI[2]),RR,paste(RR_CI[1],RR_CI[2]))

#----------------------------------------------------------
#3. Selection bias- Type II (generalizability bias)- english speaking 
#---------------------------------------------------------

set.seed(294756) #for replicability 

#----- Step1 : obtain  bias parameter estimates 

E_0mu <- log(5)
E_0sd <- (log(9)-log(4))/(2*1.96)

E_xmu <- log(0.88)
E_xsd <- (log(0.93)-log(0.84))/(2*1.96)

E_ymu <- log(1.2)
E_ysd <- (log(1.5)-log(0.90))/(2*1.96)

# S_xymu <-log(1.3)
# S_xysd <- (log(1.6)-log(1.1))/(2*1.96)

E_xymu <-log(0.70)
E_xysd <- (log(0.80)-log(0.5))/(2*1.96)

#hist(exp(rnorm(10000,S_xymu,S_xysd)))



#vector to save results 
coef_RD <- c()
coef_RR <- c()


#no of bootstrap samples
I <- 1000

#step 2: fit the bias model 

for(i in 1:I){
  
  print(i)
  bootstrap_indices <- sample(1:nrow(dat), nrow(dat), replace=TRUE)


E_0 <- rnorm(1,S_0mu,S_0sd)
E_x <- rnorm(1,S_xmu,S_xsd)
E_y <- rnorm(1,S_ymu,S_ysd)
E_xy <- rnorm(1,S_xymu,S_xysd)

df <- dat[bootstrap_indices,]

df[cat.vars]=lapply(df[cat.vars],as.numeric)
df[cat.vars]= df[cat.vars]-1

df <-  df %>% 
  mutate( pE = expit(E_0 +E_x  * any_brfeed_an + E_y * asthma_dx_current_6y + E_xy* asthma_dx_current_6y*any_brfeed_an )
  )

cat.vars_an <- names(df)[c(3:7,9:17)]
df[cat.vars_an]=lapply(df[cat.vars_an],factor)



outcome <- c("asthma_dx_current_6y")
exp_confounders <- c("any_brfeed_an","sexii","mode_delivery","pre_term","LBW_infants",
                     "mother_age_atbirth","smoke_mother_pregnancyii","parent_cob_3cat",
                     "SEIFA_quint" ,"mother_asthmaii","father_asthmaii","family_history_allergyii",
                     "number_sibs")   


df$ipw <- 1/df$pE

df <- df %>% filter(!is.na(ipw))
final_RD <- svyglm(as.formula(paste(outcome, paste(exp_confounders, collapse=" + "), sep=" ~ ")), family=gaussian("identity"),
                     design=svydesign(id=~1, weights=~ipw, data=df))

coef_RD[i] <- coef(final_RD)[2]

#-------------------Risk ratio
final_RR <- svyglm(as.formula(paste(outcome, paste(exp_confounders, collapse=" + "), sep=" ~ ")), family=poisson("log"),
                   design=svydesign(id=~1, weights=~ipw, data=df))

coef_RR[i] <- coef(final_RR)[2]

}



RD <- round(median(coef_RD),2) #estimate
RD_CI <- round(quantile(coef_RD,c(.025, .975)),2) # 95% CI


RR <- round(exp(median(coef_RR)),2) # estimate
RR_CI <- round(exp(quantile(coef_RR,c(.025, .975))),2) # 95% CI 

Table5[5,] <- c("SB-generalizability- ethnicity",RD,
                paste(RD_CI[1],RD_CI[2]),RR,paste(RR_CI[1],RR_CI[2]))

#------------------------------------------------
#4. Measurement bias (in exposure )
#------------------------------------------------

set.seed(738929) #for replicability 

#----- Step1 : obtain  bias parameter estimates 

#sensitivity
mu_sensi <- 0.84
sd_sensi <- (0.90-0.82)/(2*1.96)

n <- mu_sensi*(1-mu_sensi)
d <- (sd_sensi^2)
sensi_a <- mu_sensi*((n/d)-1)
sensi_b<-  (1-mu_sensi)*((n/d)-1)

hist(rbeta(10000,sensi_a,sensi_b))


#specificity
mu_speci <- 0.84
sd_speci <- (0.85-0.70)/(2*1.96)

n <- mu_speci*(1-mu_speci)
d <- (sd_speci^2)
speci_a <- mu_speci*((n/d)-1)
speci_b<-  (1-mu_speci)*((n/d)-1)


hist(rbeta(10000,speci_a,speci_b))

#population BF initiation prevalence
 # mu_p <- log(0.85)
 # mu_sd <- (log(0.90)-log(0.83))/(2*1.96)
 # hist(exp(rnorm(10000,mu_p,mu_sd)))
 # 
# 
# n <- mu_p*(1-mu_p)
# d <- (mu_sd^2)
# a <- mu_p*((n/d)-1)
# b<-  (1-mu_p)*((n/d)-1)


#Other bias parameters 
exp_model <- glm(as.formula(paste("any_brfeed_an", paste(c(outcome,confounders),collapse= "+")
                                  ,sep="~")),family=binomial(logit), data=dat)


#vector to save results 
coef_RD <- c()
coef_RR <- c()


#no of bootstrap samples
I <-1000
i=1

#step 2: fit the bias model 


while (length(na.omit(coef_RD))<I){
  

  print(i)
  bootstrap_indices <- sample(1:nrow(dat), nrow(dat), replace=TRUE)
  
  #sample sensitivity and specificity
  sensi <- rbeta(1, sensi_a,sensi_b)
  speci <- rbeta(1,speci_a,speci_b)
  

  #sample population prevalence
  p_bf <- runif(1,0.8,0.9)
  
  
  TP <- sensi*(nrow(dat)*p_bf)
  FN <-(nrow(dat)*p_bf)- TP
  TN <- speci*(nrow(dat)*(1-p_bf))
  FP <- (nrow(dat)*(1-p_bf))-TN
 

  x1_0     <- log((TP+FN)/(FP+TN))
  x1_x     <- log((TP/FP)/(FN/TN))
  x1_y     <-  rnorm(1, normBP(2,exp_model)[1],normBP(2,exp_model)[2])
  x1_c1    <- rnorm(1, normBP(3,exp_model)[1],normBP(3,exp_model)[2])
  x1_c2    <-rnorm(1, normBP(4,exp_model)[1],normBP(4,exp_model)[2])
  x1_c3    <-rnorm(1, normBP(5,exp_model)[1],normBP(5,exp_model)[2])
  x1_c4    <-rnorm(1, normBP(6,exp_model)[1],normBP(6,exp_model)[2])
  x1_c5    <-rnorm(1, normBP(7,exp_model)[1],normBP(7,exp_model)[2])
  x1_c6    <-rnorm(1, normBP(8,exp_model)[1],normBP(8,exp_model)[2])
  x1_c71    <-rnorm(1, normBP(9,exp_model)[1],normBP(9,exp_model)[2])
  x1_c72    <-rnorm(1, normBP(10,exp_model)[1],normBP(10,exp_model)[2])
  x1_c81    <-rnorm(1, normBP(11,exp_model)[1],normBP(11,exp_model)[2])
  x1_c82    <-rnorm(1, normBP(12,exp_model)[1],normBP(12,exp_model)[2])
  x1_c83    <-rnorm(1, normBP(13,exp_model)[1],normBP(13,exp_model)[2])
  x1_c84    <-rnorm(1, normBP(14,exp_model)[1],normBP(14,exp_model)[2])
  x1_c9    <-rnorm(1, normBP(15,exp_model)[1],normBP(15,exp_model)[2])
  x1_c10   <-rnorm(1, normBP(16,exp_model)[1],normBP(16,exp_model)[2])
  x1_c11   <-rnorm(1, normBP(17,exp_model)[1],normBP(17,exp_model)[2])
  x1_c121   <- rnorm(1, normBP(18,exp_model)[1],normBP(18,exp_model)[2])
  x1_c122   <-rnorm(1, normBP(19,exp_model)[1],normBP(19,exp_model)[2])
  x1_c123  <-rnorm(1, normBP(20,exp_model)[1],normBP(20,exp_model)[2])

 
  
  df <- dat[bootstrap_indices,]
  n=nrow(dat)
  
  #code to numeric for rbinom
  cat.vars2 <- c("any_brfeed_an","sexii","mode_delivery","pre_term")
  df[cat.vars2]=lapply(df[cat.vars2],as.numeric)
  
  df[cat.vars]=lapply(df[cat.vars],as.numeric)
  df[cat.vars]= df[cat.vars]-1
  
  

df <-  df %>% 
  mutate( x_pred = rbinom(n, 1, expit(x1_0 + x1_x * any_brfeed_an +x1_y * asthma_dx_current_6y + x1_c1* sexii +
                                        x1_c2* mode_delivery+ x1_c3*pre_term+ x1_c4*LBW_infants +
                                        x1_c5*mother_age_atbirth+x1_c6*smoke_mother_pregnancyii +
                                        x1_c71*parent_cob_3cat_1_or_both_asian +
                                        x1_c72*parent_cob_3cat_other +x1_c81*SEIFA_quint_2 +
                                        x1_c82*SEIFA_quint_3 + x1_c83*SEIFA_quint_4 +
                                        x1_c84*SEIFA_quint_5 + x1_c9*mother_asthmaii +
                                        x1_c10*father_asthmaii + x1_c11*family_history_allergyii +
                                        x1_c121*number_sibs_1_sibling + x1_c122*number_sibs_2_siblings +
                                        x1_c123*number_sibs_3_or_more_siblings )))


df$x_pred_cat <- as.factor(df$x_pred)
levels(df$x_pred_cat) <-c("0","1") 

   p<- table(df$asthma_dx_current_6y,df$x_pred_cat)
   if(p[1,1]>0 & p[1,2]>0 & p[2,1]>0 & p[2,2]>0){
  


#recode to factor
cat.vars_an <- names(df)[c(3:7,9:17)]
df[cat.vars_an]=lapply(df[cat.vars_an],factor)


confounders <- c("sexii" ,"mode_delivery","pre_term","LBW_infants",
                 "mother_age_atbirth","smoke_mother_pregnancyii","parent_cob_3cat","SEIFA_quint" ,"mother_asthmaii","father_asthmaii","family_history_allergyii",
                 "number_sibs")

final_RD <- glm(as.formula(paste("asthma_dx_current_6y", paste(c("x_pred",confounders),collapse= "+"),sep="~")),data=df,family = gaussian("identity"))




coef_RD[i] <- coef(final_RD)[2]

#-------Risk ratio

final_RR <- glm(as.formula(paste("asthma_dx_current_6y", paste(c("x_pred",confounders),collapse= "+"),sep="~")),data=df,family = poisson("log"))

coef_RR[i] <- coef(final_RR)[2]

   } else {coef_RD[i] <- NA
   coef_RD[i] <- NA}

  i <- i+1
}



RD <- round(median(coef_RD,na.rm = T),2) #estimate
RD_CI <- round(quantile(coef_RD,c(.025, .975),na.rm = T),2) # 95% CI


RR <- round(exp(median(coef_RR,na.rm = T)),2) # estimate
RR_CI <- round(exp(quantile(coef_RR,c(.025, .975),na.rm = T)),2) # 95% CI 

Table5[6,] <- c("MB-X",RD,
                paste(RD_CI[1],RD_CI[2]),RR,paste(RR_CI[1],RR_CI[2]))


#------------------------------------------------
#5. Measurement bias (in outcome )
#------------------------------------------------

set.seed(86424) #for replicability 

#----- Step1 : obtain  bias parameter estimates 


#sensitivity
mu_sensi <- 0.80
sd_sensi <- (0.90-0.77)/(2*1.96)

n <- mu_sensi*(1-mu_sensi)
d <- (sd_sensi^2)
sensi_c <- mu_sensi*((n/d)-1)
sensi_d<-  (1-mu_sensi)*((n/d)-1)

hist(rbeta(10000,sensi_c,sensi_d))


#specificity
mu_speci <- 0.90
sd_speci <- (0.95-0.81)/(2*1.96)

n <- mu_speci*(1-mu_speci)
d <- (sd_speci^2)
speci_c <- mu_speci*((n/d)-1)
speci_d<-  (1-mu_speci)*((n/d)-1)

hist(rbeta(10000,speci_c,speci_d))

#population asthma prevalence
c <- 0.21
d <- 0.30


#Other bias parameters 
y_model <-  glm(as.formula(paste("asthma_dx_current_6y", paste(c("any_brfeed_an",confounders),collapse= "+"),sep="~")),family=binomial(logit), data=dat)

#vector to save results 
coef_RD <- c()
coef_RR <- c()


#no of bootstrap samples
I <- 1000

#step 2: fit the bias model 

for(i in 1:I){
  
  print(i)
  bootstrap_indices <- sample(1:nrow(dat), nrow(dat), replace=TRUE)
  
  df <-  dat[bootstrap_indices,]
  
  #sample sensitivity and specificity
  sensi <- rbeta(1, sensi_c,sensi_d)
  speci <- rbeta(1,speci_c,speci_d)
  
  
  #sample population prevalence
  p_asth <-runif(1,c,d)
  
  TP <- sensi*(nrow(dat)*p_asth)
  FN <-(nrow(dat)*p_asth)- TP
  TN <- speci*(nrow(dat)*(1-p_asth))
  FP <- (nrow(dat)*(1-p_asth))-TN
  
  
  y1_0     <- log((TP+FN)/(FP+TN))
  y1_y     <- log((TP/FP)/(FN/TN))
y1_xpred <-  rnorm(1, normBP(2,y_model)[1],normBP(2,y_model)[2])
y1_c1    <- rnorm(1, normBP(3,y_model)[1],normBP(3,y_model)[2])
y1_c2    <- rnorm(1, normBP(4,y_model)[1],normBP(4,y_model)[2])
y1_c3    <- rnorm(1, normBP(5,y_model)[1],normBP(5,y_model)[2])
y1_c4    <- rnorm(1, normBP(6,y_model)[1],normBP(6,y_model)[2])
y1_c5    <- rnorm(1, normBP(7,y_model)[1],normBP(7,y_model)[2])
y1_c6    <- rnorm(1, normBP(8,y_model)[1],normBP(8,y_model)[2])
y1_c71    <-rnorm(1, normBP(9,y_model)[1],normBP(9,y_model)[2])
y1_c72    <- rnorm(1, normBP(10,y_model)[1],normBP(10,y_model)[2])
y1_c81    <- rnorm(1, normBP(11,y_model)[1],normBP(11,y_model)[2])
y1_c82    <- rnorm(1, normBP(12,y_model)[1],normBP(12,y_model)[2])
y1_c83    <- rnorm(1, normBP(13,y_model)[1],normBP(13,y_model)[2])
y1_c84    <- rnorm(1, normBP(14,y_model)[1],normBP(14,y_model)[2])
y1_c9    <- rnorm(1, normBP(15,y_model)[1],normBP(15,y_model)[2])
y1_c10   <- rnorm(1, normBP(16,y_model)[1],normBP(16,y_model)[2])
y1_c11   <- rnorm(1, normBP(17,y_model)[1],normBP(17,y_model)[2])
y1_c121   <-rnorm(1, normBP(18,y_model)[1],normBP(18,y_model)[2])
y1_c122   <- rnorm(1, normBP(19,y_model)[1],normBP(19,y_model)[2])
y1_c123   <- rnorm(1, normBP(20,y_model)[1],normBP(20,y_model)[2])




#code to numeric for rbinom
cat.vars2 <- c("any_brfeed_an","sexii","mode_delivery","pre_term")
df[cat.vars2]=lapply(df[cat.vars2],as.numeric)

df[cat.vars]=lapply(df[cat.vars],as.numeric)
df[cat.vars]= df[cat.vars]-1


n=nrow(dat)

df <- df %>% 
  mutate(y_pred=rbinom(n, 1, expit(y1_0 +y1_y * asthma_dx_current_6y +  y1_xpred*any_brfeed_an+
                                     y1_c1 * sexii+y1_c2 * mode_delivery+
                                y1_c3 * pre_term+y1_c4 * LBW_infants+
                                  y1_c5 * mother_age_atbirth+ y1_c6 * smoke_mother_pregnancyii+
                                  y1_c71 *parent_cob_3cat_1_or_both_asian +
                                y1_c71*parent_cob_3cat_other+y1_c81*SEIFA_quint_2
                                +y1_c82*SEIFA_quint_3+ y1_c83*SEIFA_quint_4+
                                  y1_c84*SEIFA_quint_5+y1_c9* mother_asthmaii+
                                y1_c10 * father_asthmaii+y1_c11 * family_history_allergyii+
                                  y1_c121*number_sibs_1_sibling+y1_c122*number_sibs_2_siblings+
                                y1_c123*number_sibs_3_or_more_siblings ))
  )



#recode to factor
cat.vars_an <- names(df)[c(3:7,9:17)]
df[cat.vars_an]=lapply(df[cat.vars_an],factor)

confounders <- c("sexii" ,"mode_delivery","pre_term","LBW_infants",
                 "mother_age_atbirth","smoke_mother_pregnancyii","parent_cob_3cat","SEIFA_quint" ,"mother_asthmaii","father_asthmaii","family_history_allergyii",
                 "number_sibs")


final_RD <- glm(as.formula(paste("y_pred", paste(c("any_brfeed_an",confounders),collapse= "+"),sep="~")),data=df,family = gaussian("identity"))
coef_RD[i] <- coef(final_RD)[2]

#-------Risk ratio

final_RR <- glm(as.formula(paste("y_pred", paste(c("any_brfeed_an",confounders),collapse= "+"),sep="~")),data=df,family = poisson("log"))

coef_RR[i] <- coef(final_RR)[2]

}


RD <- round(median(coef_RD),2) #estimate
RD_CI <- round(quantile(coef_RD,c(.025, .975)),2) # 95% CI


RR <- round(exp(median(coef_RR)),2) # estimate
RR_CI <- round(exp(quantile(coef_RR,c(.025, .975))),2) # 95% CI 

Table5[7,] <- c("MB-Y",RD,
                paste(RD_CI[1],RD_CI[2]),RR,paste(RR_CI[1],RR_CI[2]))


#-----------------------------
#6. Confounding bias 
#----------------------------

set.seed(597205) #for replicability 

#----- Step1 : obtain  bias parameter estimates 

mu_u0 <- log(0.12)
sd_u0 <-  (log(0.20)-log(0.07))/(2*1.96)


mu_ux <- log(0.7)
sd_ux <-  (log(0.8)-log(0.5))/(2*1.96)

mu_uy <- log(1.1)
sd_uy <-  (log(1.2)-log(0.88))/(2*1.96)


#vector to save results 
coef_RD <- c()
coef_RR <- c()


#no of bootstrap samples
I <- 1000

#step 2: fit the bias model 

for(i in 1:I){
  
  print(i)
  bootstrap_indices <- sample(1:nrow(dat), nrow(dat), replace=TRUE)
  
  u1_0     <- rnorm(1,mean=mu_u0,sd=sd_u0)
  u1_x     <- rnorm(1,mean=mu_ux,sd=sd_ux)
  u1_y     <- rnorm(1,mean=mu_uy,sd=sd_uy)
  
  df <- dat[bootstrap_indices,]  
    
  
  #code to numeric for rbinom
  cat.vars2 <- c("any_brfeed_an","sexii","mode_delivery","pre_term")
  df[cat.vars2]=lapply(df[cat.vars2],as.numeric)
  
  df[cat.vars]=lapply(df[cat.vars],as.numeric)
  df[cat.vars]= df[cat.vars]-1
  
  df <- df %>% 
    mutate(u1pred = rbinom(nrow(dat),1, expit(u1_0+u1_x*any_brfeed_an+u1_y*asthma_dx_current_6y)))
  
  
  #recode to factor
  cat.vars_an <- names(df)[c(3:7,9:17)]
  df[cat.vars_an]=lapply(df[cat.vars_an],factor)
    df$u1pred <- as.factor(df$u1pred)


exp_confounders <- c("any_brfeed_an","sexii","mode_delivery","pre_term","LBW_infants",
                     "mother_age_atbirth","smoke_mother_pregnancyii","parent_cob_3cat",
                     "SEIFA_quint" ,"mother_asthmaii","father_asthmaii","family_history_allergyii",
                     "number_sibs","u1pred")  


outcome <- c("asthma_dx_current_6y")

final_RD <- glm(as.formula(paste(outcome, paste(exp_confounders, collapse=" + "), sep=" ~ ")),
            family = gaussian("identity"), df)

coef_RD[i] <- coef(final_RD)[2]
#---------

final_RR <- glm(as.formula(paste(outcome, paste(exp_confounders, collapse=" + "), sep=" ~ ")),
            family = poisson("log"), df)

coef_RR[i] <- coef(final_RR)[2]
}



RD <- round(median(coef_RD),2) #estimate
RD_CI <- round(quantile(coef_RD,c(.025, .975)),2) # 95% CI


RR <- round(exp(median(coef_RR)),2) # estimate
RR_CI <- round(exp(quantile(coef_RR,c(.025, .975))),2) # 95% CI 

Table5[8,] <- c("CB",RD,
                paste(RD_CI[1],RD_CI[2]),RR,paste(RR_CI[1],RR_CI[2]))

#-----------------------------------------------------------------
# 7. Adjusting for all biases simultaneously using the imputation approach 
#------------------------------------------------------------------

set.seed(3920124) #for replicability 

#vector to save results 
coef_RD <- c()
coef_RR <- c()

#----- Step1 : obtain  bias parameter estimates 
mu_xpred <-log(0.6) 
sd_xpred <- (log(0.8)-log(0.5))/(2*1.96)

#hist(exp(rnorm(100000,mu_xpred,sd_xpred)))
  
#all other parameters specified above


#no of bootstrap samples
I <- 1000
i=1

#step 2: fit the bias model 

while (length(na.omit(coef_RD))<I){
  
  print(i)
  bootstrap_indices <- sample(1:nrow(dat), nrow(dat), replace=TRUE)
  df <- dat[bootstrap_indices,]  
  
  
  #selection bias 1
  s1_0     <- rnorm(1,normBP(1,propmodel)[1],normBP(1,propmodel)[2])
  s1_x     <- rnorm(1, normBP(2,propmodel)[1],normBP(2,propmodel)[2])
  
  s1_y     <- rnorm(1,mean=mu_y,sd=sd_y)
  
  s1_c1    <-rnorm(1, normBP(3,propmodel)[1],normBP(3,propmodel)[2])
  s1_c2    <-rnorm(1, normBP(4,propmodel)[1],normBP(4,propmodel)[2])
  s1_c3    <-rnorm(1, normBP(5,propmodel)[1],normBP(5,propmodel)[2])
  s1_c4    <-rnorm(1, normBP(6,propmodel)[1],normBP(6,propmodel)[2])
  s1_c5    <-rnorm(1, normBP(7,propmodel)[1],normBP(7,propmodel)[2])
  s1_c6    <-rnorm(1, normBP(8,propmodel)[1],normBP(8,propmodel)[2])
  s1_c71    <-rnorm(1, normBP(9,propmodel)[1],normBP(9,propmodel)[2])
  s1_c72    <-rnorm(1, normBP(10,propmodel)[1],normBP(10,propmodel)[2])
  
  s1_c81    <-rnorm(1, normBP(11,propmodel)[1],normBP(11,propmodel)[2])
  s1_c82    <-rnorm(1, normBP(12,propmodel)[1],normBP(12,propmodel)[2])
  s1_c83    <-rnorm(1, normBP(13,propmodel)[1],normBP(13,propmodel)[2])
  s1_c84    <-rnorm(1, normBP(14,propmodel)[1],normBP(14,propmodel)[2])
  
  s1_c9    <-rnorm(1, normBP(15,propmodel)[1],normBP(15,propmodel)[2])
  s1_c10   <-rnorm(1, normBP(16,propmodel)[1],normBP(16,propmodel)[2])
  s1_c11   <-rnorm(1, normBP(17,propmodel)[1],normBP(17,propmodel)[2])
  
  s1_c121   <-rnorm(1, normBP(18,propmodel)[1],normBP(18,propmodel)[2])
  s1_c122   <-rnorm(1, normBP(19,propmodel)[1],normBP(19,propmodel)[2])
  s1_c123   <-rnorm(1, normBP(20,propmodel)[1],normBP(20,propmodel)[2])
  
  
  #selection bias 2 (consent and ethnicity)
  S_0 <- rnorm(1,S_0mu,S_0sd)
  S_x <- rnorm(1,S_xmu,S_xsd)
  S_y <- rnorm(1,S_ymu,S_ysd)
  S_xy <- rnorm(1,S_xymu,S_xysd)
  
  
  E_0 <- rnorm(1,S_0mu,S_0sd)
  E_x <- rnorm(1,S_xmu,S_xsd)
  E_y <- rnorm(1,S_ymu,S_ysd)
  E_xy <- rnorm(1,S_xymu,S_xysd)
  
  
  
  #sample sensitivity and specificity
  sensi <- rbeta(1, sensi_a,sensi_b)
  speci <- rbeta(1,speci_a,speci_b)
  
  
  #sample population prevalence
  # p_bf <- rbeta(1, a,b)
  p_bf <- runif(1,0.8,0.90)
  
  
  TP <- sensi*(nrow(dat)*p_bf)
  FN <-(nrow(dat)*p_bf)- TP
  TN <- speci*(nrow(dat)*(1-p_bf))
  FP <- (nrow(dat)*(1-p_bf))-TN
  
  #exposure misclassification
  x1_0     <- log((TP+FN)/(FP+TN))
  x1_x     <- log((TP/FP)/(FN/TN))
  x1_y     <-  rnorm(1, normBP(2,exp_model)[1],normBP(2,exp_model)[2])
  x1_c1    <- rnorm(1, normBP(3,exp_model)[1],normBP(3,exp_model)[2])
  x1_c2    <-rnorm(1, normBP(4,exp_model)[1],normBP(4,exp_model)[2])
  x1_c3    <-rnorm(1, normBP(5,exp_model)[1],normBP(5,exp_model)[2])
  x1_c4    <-rnorm(1, normBP(6,exp_model)[1],normBP(6,exp_model)[2])
  x1_c5    <-rnorm(1, normBP(7,exp_model)[1],normBP(7,exp_model)[2])
  x1_c6    <-rnorm(1, normBP(8,exp_model)[1],normBP(8,exp_model)[2])
  x1_c71    <-rnorm(1, normBP(9,exp_model)[1],normBP(9,exp_model)[2])
  x1_c72    <-rnorm(1, normBP(10,exp_model)[1],normBP(10,exp_model)[2])
  x1_c81    <-rnorm(1, normBP(11,exp_model)[1],normBP(11,exp_model)[2])
  x1_c82    <-rnorm(1, normBP(12,exp_model)[1],normBP(12,exp_model)[2])
  x1_c83    <-rnorm(1, normBP(13,exp_model)[1],normBP(13,exp_model)[2])
  x1_c84    <-rnorm(1, normBP(14,exp_model)[1],normBP(14,exp_model)[2])
  x1_c9    <-rnorm(1, normBP(15,exp_model)[1],normBP(15,exp_model)[2])
  x1_c10   <-rnorm(1, normBP(16,exp_model)[1],normBP(16,exp_model)[2])
  x1_c11   <-rnorm(1, normBP(17,exp_model)[1],normBP(17,exp_model)[2])
  x1_c121   <- rnorm(1, normBP(18,exp_model)[1],normBP(18,exp_model)[2])
  x1_c122   <-rnorm(1, normBP(19,exp_model)[1],normBP(19,exp_model)[2])
  x1_c123  <-rnorm(1, normBP(20,exp_model)[1],normBP(20,exp_model)[2])

  
  #outcome misclassification 
  #sample sensitivity and specificity
  sensi <- rbeta(1, sensi_c,sensi_d)
  speci <- rbeta(1,speci_c,speci_d)
  
  
  #sample population prevalence
  p_asth <-runif(1,c,d)
  
  TP <- sensi*(nrow(dat)*p_asth)
  FN <-(nrow(dat)*p_asth)- TP
  TN <- speci*(nrow(dat)*(1-p_asth))
  FP <- (nrow(dat)*(1-p_asth))-TN
  
  
  y1_0     <- log((TP+FN)/(FP+TN))
  y1_y     <- log((TP/FP)/(FN/TN))
  y1_xpred <-  rnorm(1, normBP(2,y_model)[1],normBP(2,y_model)[2])
  y1_c1    <- rnorm(1, normBP(3,y_model)[1],normBP(3,y_model)[2])
  y1_c2    <- rnorm(1, normBP(4,y_model)[1],normBP(4,y_model)[2])
  y1_c3    <- rnorm(1, normBP(5,y_model)[1],normBP(5,y_model)[2])
  y1_c4    <- rnorm(1, normBP(6,y_model)[1],normBP(6,y_model)[2])
  y1_c5    <- rnorm(1, normBP(7,y_model)[1],normBP(7,y_model)[2])
  y1_c6    <- rnorm(1, normBP(8,y_model)[1],normBP(8,y_model)[2])
  y1_c71    <-rnorm(1, normBP(9,y_model)[1],normBP(9,y_model)[2])
  y1_c72    <- rnorm(1, normBP(10,y_model)[1],normBP(10,y_model)[2])
  y1_c81    <- rnorm(1, normBP(11,y_model)[1],normBP(11,y_model)[2])
  y1_c82    <- rnorm(1, normBP(12,y_model)[1],normBP(12,y_model)[2])
  y1_c83    <- rnorm(1, normBP(13,y_model)[1],normBP(13,y_model)[2])
  y1_c84    <- rnorm(1, normBP(14,y_model)[1],normBP(14,y_model)[2])
  y1_c9    <- rnorm(1, normBP(15,y_model)[1],normBP(15,y_model)[2])
  y1_c10   <- rnorm(1, normBP(16,y_model)[1],normBP(16,y_model)[2])
  y1_c11   <- rnorm(1, normBP(17,y_model)[1],normBP(17,y_model)[2])
  y1_c121   <-rnorm(1, normBP(18,y_model)[1],normBP(18,y_model)[2])
  y1_c122   <- rnorm(1, normBP(19,y_model)[1],normBP(19,y_model)[2])
  y1_c123   <- rnorm(1, normBP(20,y_model)[1],normBP(20,y_model)[2])

  y1_xpred <- rnorm(1,mu_xpred,sd_xpred)
  
  
  #confounding 
  u1_0     <- rnorm(1,mean=mu_u0,sd=sd_u0)
  u1_x     <- rnorm(1,mean=mu_ux,sd=sd_ux)
  u1_y     <- rnorm(1,mean=mu_uy,sd=sd_uy)
  
  
  #code to numeric for rbinom
  cat.vars2 <- c("any_brfeed_an","sexii","mode_delivery","pre_term")
  df[cat.vars2]=lapply(df[cat.vars2],as.numeric)
  
  df[cat.vars]=lapply(df[cat.vars],as.numeric)
  df[cat.vars]= df[cat.vars]-1

  
  n=nrow(df)
  
df <- df %>% 
  mutate( x_pred = rbinom(n, 1, expit(x1_0 + x1_x * any_brfeed_an+ x1_c1* sexii+x1_c2* mode_delivery+x1_c3+pre_term+x1_c4+LBW_infants+
                                       x1_c5*mother_age_atbirth+x1_c6*smoke_mother_pregnancyii+x1_c71*parent_cob_3cat_1_or_both_asian+x1_c72*parent_cob_3cat_other+
                                       x1_c81*SEIFA_quint_2+x1_c82*SEIFA_quint_3+x1_c83*SEIFA_quint_4+x1_c84*SEIFA_quint_5+x1_c9*mother_asthmaii+
                                       x1_c10*father_asthmaii+x1_c11*family_history_allergyii+x1_c121*number_sibs_1_sibling+x1_c122*number_sibs_2_siblings+
                                       x1_c123*number_sibs_3_or_more_siblings)),
         y_pred=rbinom(n, 1, expit(y1_0 +y1_y * asthma_dx_current_6y +   y1_c1 * sexii+y1_c2 * mode_delivery+
                                     y1_c3 * pre_term+y1_c4 * LBW_infants+y1_c5 * mother_age_atbirth+
                                     y1_c6 * smoke_mother_pregnancyii+y1_c71 * parent_cob_3cat_1_or_both_asian+y1_c71*parent_cob_3cat_other+
                                     y1_c81*SEIFA_quint_2+y1_c82*SEIFA_quint_3+y1_c83*SEIFA_quint_4+y1_c84*SEIFA_quint_5+y1_c9* mother_asthmaii+
                                     y1_c10 * father_asthmaii+
                                     y1_c11 * family_history_allergyii+y1_c121*number_sibs_1_sibling+y1_c122*number_sibs_2_siblings+
                                     y1_c123*number_sibs_3_or_more_siblings+
                                   y1_xpred*x_pred )),
         u1pred = rbinom(nrow(dat),1, expit(u1_0+u1_x*x_pred+u1_y*y_pred)),
         pS1 = expit(s1_0 + s1_x * any_brfeed_an + s1_y * asthma_dx_current_6y + s1_c1 * sexii+s1_c2 * mode_delivery+
                      s1_c3 * pre_term+s1_c4 * LBW_infants+s1_c5 * mother_age_atbirth+ 
                       s1_c6 * smoke_mother_pregnancyii+s1_c71 * parent_cob_3cat_1_or_both_asian +
                      s1_c71*parent_cob_3cat_other+
                      s1_c81* SEIFA_quint_2+s1_c82* SEIFA_quint_3+ s1_c83* SEIFA_quint_4+s1_c84* SEIFA_quint_5 +
                       s1_c9* mother_asthmaii+s1_c10 * father_asthmaii+s1_c11 * family_history_allergyii+
                      s1_c121 * number_sibs_1_sibling +s1_c122 *number_sibs_2_siblings +s1_c122* number_sibs_3_or_more_siblings),
         pS2 = expit(S_0 +S_x  * any_brfeed_an + S_y * asthma_dx_current_6y + S_xy* asthma_dx_current_6y*any_brfeed_an ),
         pS3 = expit(E_0 +E_x  * any_brfeed_an + E_y * asthma_dx_current_6y + E_xy* asthma_dx_current_6y*any_brfeed_an )
         
       
  )




df$x_pred_cat <- as.factor(df$x_pred)
levels(df$x_pred_cat) <-c("0","1")

p<- table(df$asthma_dx_current_6y,df$x_pred_cat)
if(p[1,1]>0 & p[1,2]>0 & p[2,1]>0 & p[2,2]>0){


df$PS1S2S3 <- df$pS1*df$pS2*df$pS3

#recode to factor
cat.vars_an <- names(df)[c(3:7,9:17)]
df[cat.vars_an]=lapply(df[cat.vars_an],factor)
df$u1pred <- as.factor(df$u1pred)



df$ipw <- 1/df$PS1S2S3

df <- df %>% filter(!is.na(ipw))

final_RD<- svyglm(as.formula(paste("y_pred",paste(c("x_pred",confounders,"u1pred"), collapse=" + "), sep=" ~ ")), family=gaussian("identity"),
                     design=svydesign(id=~1, weights=~ipw, data=df))
coef_RD[i] <- coef(final_RD)[2]

final_RR<- svyglm(as.formula(paste("y_pred",paste(c("x_pred",confounders,"u1pred"), collapse=" + "), sep=" ~ ")), family=poisson("log"),
                design=svydesign(id=~1, weights=~ipw, data=df))

coef_RR[i] <- coef(final_RR)[2]


} else {coef_RD[i] <- NA
coef_RR[i] <- NA}

i <- i+1

}



RD <- round(median(coef_RD,na.rm = T),2) #estimate
RD_CI <- round(quantile(coef_RD,c(.025, .975),na.rm = T),2) # 95% CI


RR <- round(exp(median(coef_RR,na.rm = T)),2) # estimate
RR_CI <- round(exp(quantile(coef_RR,c(.025, .975),na.rm = T)),2) # 95% CI 
Table5[2,] <- c("All biases",RD,
                paste(RD_CI[1],RD_CI[2]),RR,paste(RR_CI[1],RR_CI[2]))


#order the table


write.csv (Table5,"case study results.csv")
#-------------------------------Figure 6

Table_RD <- data.frame(matrix(ncol = 7, nrow = 8))
names(Table_RD) <- c("Biases adjusted for","Risk difference","LCIRD","UCIRD","Risk ratio","LCIRR","UCIRR")
Table_RD$`Biases adjusted for` <- c("Primary analysis","SB-collider stratification","SB-generalizability (consent)",
                                    "SB-generalizability (ethnicity)","MB-A","MB-Y","CB","All-biases")
Table_RD$`Risk difference` <- c(Table5$`Risk difference`)
Table_RD$`Risk difference` <- as.numeric(Table_RD$`Risk difference` )


Table_RD$`Risk ratio` <- c(Table5$`Risk ratio`)
Table_RD$`Risk ratio` <- as.numeric(Table_RD$`Risk ratio` )


Table_RD$UCIRD <- c(0.02,0.03,0.03,0.03,0.07,0.06,0.04,0.05)

Table_RD$LCIRD <- c((-0.11),(-0.09),(-0.10),(-0.10),(-0.19),(-0.19),(-0.10),(-0.26))

Table_RD$UCIRR <- c(1.15,1.37,1.31,1.29,1.92,1.14,1.37,1.14)

Table_RD$LCIRR<- c(0.82,0.56,0.54,0.53,0.41,0.64,0.57,0.54)

str(Table_RD)

Table_RD$`Biases adjusted for` <- factor(Table_RD$`Biases adjusted for`,
                                         levels=c("Primary analysis","All-biases","SB-generalizability (consent)","SB-generalizability (ethnicity)",
                                                  "SB-collider stratification","MB-A","MB-Y","CB"
                                         ))

ggplot(Table_RD, aes(x=`Biases adjusted for`, y=`Risk difference`,color=`Biases adjusted for` )) + geom_pointrange(aes(ymin = LCIRD, ymax = UCIRD))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=13))


ggplot(Table_RD, aes(x=`Biases adjusted for`, y=`Risk ratio`,color=`Biases adjusted for` )) + geom_pointrange(aes(ymin = LCIRR, ymax = UCIRR))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=13))
