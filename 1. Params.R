##################################################################################
#  Obtain parameters for the data generation 
#  Data used: HealthNuts data 
#  Author: Rushani Wijesuriya 22nd of July 2023                             #
##################################################################################


rm (list=ls())

setwd("C:/Users/rushani.wijesuriya/OneDrive - Murdoch Children's Research Institute/QBA-methods paper/Simulation study")

library(readstata13)
library(dplyr)



#import data 
dat <- read.dta13("HN_BF_asthma_clean.dta", nonint.factors = TRUE,generate.factors = T)
dat <- dat %>% select("id" ,"asthma_dx_current_6y" ,"sexii" ,"mode_delivery","pre_term","birth_weight",
                      "any_brfeed_an","mother_age_atbirth","number_sibs","smokingii","smoke_mother_pregnancyii",
                      "SEIFA_quint" ,"parent_cob_3cat","mother_asthmaii","father_asthmaii","family_history_allergyii",
                      "ccare","quest_type_h3","brfeed")

#generate birth weight variable 
dat$LBW_infants <- ifelse(dat$birth_weight>=2500,0,1)

#change all factor variables to correct type
cat.vars <- c("asthma_dx_current_6y","smokingii","smoke_mother_pregnancyii","SEIFA_quint","mother_asthmaii","father_asthmaii",
              "family_history_allergyii","ccare","LBW_infants" )

dat[cat.vars]=lapply(dat[cat.vars],factor)

dat$any_brfeed_an <- relevel(dat$any_brfeed_an,ref="no")


#-------------------------------------
#summarize the observed baseline confounders: family history of asthma, family history of allergic disease,
#parental ethnicity, number of older siblings, childcare attendance, child’s sex,
#and SES. 


#1. Father asthma
F_A_no <- prop.table(table(dat$father_asthmaii))[1]
F_A_yes <- prop.table(table(dat$father_asthmaii))[2]

#2. Mother asthma
M_A_no <- prop.table(table(dat$mother_asthmaii))[1]
M_A_yes <- prop.table(table(dat$mother_asthmaii))[2]

#3. Family history of allergic disease
Fx_al_no <- prop.table(table(dat$family_history_allergyii))[1]
Fx_al_yes <- prop.table(table(dat$family_history_allergyii))[2]

#4. Parental ethnicity
CoB_aus <- prop.table(table(dat$parent_cob_3cat))[1]
CoB_asia <- prop.table(table(dat$parent_cob_3cat))[2]
CoB_other <- prop.table(table(dat$parent_cob_3cat))[3]

#5. number of older siblings
sib_no <- prop.table(table(dat$number_sibs))[1]
sib_1 <- prop.table(table(dat$number_sibs))[2]
sib_2 <- prop.table(table(dat$number_sibs))[3]
sib_3 <- prop.table(table(dat$number_sibs))[4]

#6. Child's sex
sex_m <- prop.table(table(dat$sexii))[1]
sex_f <- prop.table(table(dat$sexii))[2]

#7. SES
SES_1 <- prop.table(table(dat$SEIFA_quint))[1]
SES_2 <- prop.table(table(dat$SEIFA_quint))[2]
SES_3 <- prop.table(table(dat$SEIFA_quint))[3]
SES_4 <- prop.table(table(dat$SEIFA_quint))[4]
SES_5 <- prop.table(table(dat$SEIFA_quint))[5]


#store  parameters

vals <-  as.numeric(c(F_A_no,F_A_yes, M_A_no,M_A_yes,Fx_al_no,Fx_al_yes,CoB_aus,CoB_asia,CoB_other,sib_no,sib_1,sib_2,sib_3,
          sex_m,sex_f, SES_1,SES_2,SES_3,SES_4,SES_5))
para_name <-  c("F_A_no","F_A_yes","M_A_no","M_A_yes","Fx_al_no","Fx_al_yes","CoB_aus","CoB_asia","CoB_other",
                "sib_no","sib_1","sib_2","sib_3",
                "sex_m","sex_f", "SES_1","SES_2","SES_3","SES_4","SES_5")


#-------------------------------------------------
#maternal age at childbirth and maternal smoking 

fit <- lm(mother_age_atbirth~as.factor(SEIFA_quint),data=dat)
mage_b0 <- coefficients(fit)[1]
mage_b1 <- coefficients(fit)[2]
mage_b2 <- coefficients(fit)[3]
mage_b3 <- coefficients(fit)[4]
mage_b4 <- coefficients(fit)[5]

fit <- glm(smoke_mother_pregnancyii~as.factor(SEIFA_quint),family=binomial(link='logit'),data=dat)
msmok_b0 <- coefficients(fit)[1]
msmok_b1 <- coefficients(fit)[2]
msmok_b2 <- coefficients(fit)[3]
msmok_b3 <- coefficients(fit)[4]
msmok_b4 <- coefficients(fit)[5]

#store parameters

vals <-  c(vals,as.numeric(c(mage_b0,mage_b1, mage_b2,mage_b3,mage_b4,
                             msmok_b0,msmok_b1,msmok_b2,msmok_b3,msmok_b4)))
para_name <-c(para_name,c("mage_b0","mage_b1","mage_b2","mage_b3","mage_b4",
                          "msmok_b0","msmok_b1","msmok_b2","msmok_b3",
                "msmok_b4"))

#-------------------------------------------------
#Gestational age ~maternal age at childbirth, SES, and maternal smoking

fit <- glm(pre_term~mother_age_atbirth+as.factor(SEIFA_quint)+smoke_mother_pregnancyii,family=binomial(link='logit'),data=dat)
gest_b0 <- coefficients(fit)[1]
gest.mage_b1 <- coefficients(fit)[2]
gest.SES2_b2 <- coefficients(fit)[3]
gest.SES3_b3 <- coefficients(fit)[4]
gest.SES4_b4 <- coefficients(fit)[5]
gest.SES5_b5 <- coefficients(fit)[6]
gest.smok_b6 <- coefficients(fit)[7]


vals <-  c(vals,as.numeric(c(gest_b0,gest.mage_b1, gest.SES2_b2,gest.SES3_b3,gest.SES4_b4,gest.SES5_b5,gest.smok_b6)))
para_name <-c(para_name,c("gest_b0","gest.mage_b1","gest.SES2_b2","gest.SES3_b3","gest.SES4_b4","gest.SES5_b5","gest.smok_b6"))

#---------------------------------------------------
#mode of delivery~ on SES and maternal age at childbirth
fit <- glm(mode_delivery~mother_age_atbirth+as.factor(SEIFA_quint),family=binomial(link='logit'),data=dat)

dmode_b0 <- coefficients(fit)[1]
dmode.mage_b1 <- coefficients(fit)[2]
dmode.SES2_b2 <- coefficients(fit)[3]
dmode.SES3_b3 <- coefficients(fit)[4]
dmode.SES4_b4 <- coefficients(fit)[5]
dmode.SES5_b5 <- coefficients(fit)[6]


vals <-  c(vals,as.numeric(c(dmode_b0,dmode.mage_b1, dmode.SES2_b2,dmode.SES3_b3,dmode.SES4_b4,dmode.SES5_b5)))
para_name <-c(para_name,c("dmode_b0","dmode.mage_b1","dmode.SES2_b2","dmode.SES3_b3","dmode.SES4_b4","dmode.SES5_b5"))



#-----------------------------------------------------
#low birth weight  ~mode of delivery, gestational age, SES, maternal smoking, maternal age at childbirth 

fit <- glm(LBW_infants~mode_delivery+pre_term+mother_age_atbirth+as.factor(SEIFA_quint)+
             smokingii+smoke_mother_pregnancyii,family=binomial(link='logit'),data=dat)

LBW_b0 <- coefficients(fit)[1]
LBW.dmode_b1 <- coefficients(fit)[2]
LBW.gest_b2 <- coefficients(fit)[3]
LBW.mage_b3 <- coefficients(fit)[4]
LBW.SES2_b4 <- coefficients(fit)[5]
LBW.SES3_b5 <- coefficients(fit)[6]
LBW.SES4_b6 <- coefficients(fit)[7]
LBW.SES5_b7 <- coefficients(fit)[8]
LBW.smok_b8 <- coefficients(fit)[9]
LBW.msmok_b9 <- coefficients(fit)[10]
                  


vals <-  c(vals,as.numeric(c(LBW_b0,LBW.dmode_b1,LBW.gest_b2,LBW.mage_b3,LBW.SES2_b4,LBW.SES3_b5,
                             LBW.SES4_b6,LBW.SES5_b7,LBW.smok_b8,LBW.msmok_b9)))
para_name <-c(para_name,c("LBW_b0","LBW.dmode_b1","LBW.gest_b2","LBW.mage_b3","LBW.SES2_b4","LBW.SES3_b5",
                          "LBW.SES4_b6","LBW.SES5_b7","LBW.smok_b8","LBW.msmok_b9"))


#-----------------------------------------------------
#Exposure~mode of delivery, gestational age, SES, maternal smoking, maternal age at childbirth+
#family history of asthma, family history of allergic disease,LBW
#parental ethnicity, number of older siblings,child’s sex,


fit <- glm(any_brfeed_an~sexii+number_sibs+LBW_infants+parent_cob_3cat+mother_asthmaii+father_asthmaii+
             family_history_allergyii+
             mode_delivery+pre_term+mother_age_atbirth+as.factor(SEIFA_quint)+
             smoke_mother_pregnancyii,family=binomial(link='logit'),data=dat)


exp_b0 <- coefficients(fit)[1]
exp.sex_b1 <- coefficients(fit)[2]
exp.sibs1_b2 <- coefficients(fit)[3]
exp.sibs2_b3 <- coefficients(fit)[4]
exp.sibs3_b4 <- coefficients(fit)[5]
exp.LBW_b5 <- coefficients(fit)[6]
exp.cobasi_b6 <- coefficients(fit)[7]
exp.cobotr_b7 <- coefficients(fit)[8]
exp.msath_b8 <- coefficients(fit)[9]
exp.fsath_b9 <- coefficients(fit)[10]
exp.fhalle_b10 <- coefficients(fit)[11]
exp.dmode_b11 <- coefficients(fit)[13]
exp.gest_b12 <- coefficients(fit)[14]
exp.mage_b13 <- coefficients(fit)[15]
exp.SES2_b14<- coefficients(fit)[16]
exp.SES3_b15 <- coefficients(fit)[17]
exp.SES4_b16 <- coefficients(fit)[18]
exp.SES5_b17 <- coefficients(fit)[19]
exp.msmok_b18 <- coefficients(fit)[21]


vals <-  c(vals,as.numeric(c(exp_b0,exp.sex_b1,exp.sibs1_b2,exp.sibs2_b3,exp.sibs3_b4,exp.LBW_b5,
                             exp.cobasi_b6,exp.cobotr_b7, exp.msath_b8,exp.fsath_b9,exp.fhalle_b10,
                             exp.dmode_b11,exp.gest_b12, exp.mage_b13,exp.SES2_b14,
                             exp.SES3_b15,exp.SES4_b16,exp.SES5_b17,exp.msmok_b18)))
para_name <-c(para_name,c("exp_b0","exp.sex_b1","exp.sibs1_b2","exp.sibs2_b3","exp.sibs3_b4","exp.LBW_b5",
                          "exp.cobasi_b6","exp.cobotr_b7", "exp.msath_b8","exp.fsath_b9","exp.fhalle_b10",
                          "exp.dmode_b11","exp.gest_b12", "exp.mage_b13","exp.SES2_b14",
                          "exp.SES3_b15","exp.SES4_b16","exp.SES5_b17","exp.msmok_b18"))


param <- data.frame(para_name,vals)

#-----------------------------------------------------
#Outcome~breast feeding +mode of delivery, gestational age, SES, maternal smoking, maternal age at childbirth+
#family history of asthma, family history of allergic disease,LBW
#parental ethnicity, number of older siblings, child’s sex,

fit <- glm(asthma_dx_current_6y~any_brfeed_an+sexii+number_sibs+LBW_infants+parent_cob_3cat+mother_asthmaii+father_asthmaii+
             family_history_allergyii+
             mode_delivery+pre_term+mother_age_atbirth+as.factor(SEIFA_quint)+
             smoke_mother_pregnancyii,family=binomial(link='logit'),data=dat)

out_b0 <- coefficients(fit)[1]
out.exp_b1 <- coefficients(fit)[2]
out.sex_b2 <- coefficients(fit)[3]
out.sibs1_b3 <- coefficients(fit)[4]
out.sibs2_b4 <- coefficients(fit)[5]
out.sibs3_b5 <- coefficients(fit)[6]
out.LBW_b6 <- coefficients(fit)[7]
out.cobasi_b7 <- coefficients(fit)[8]
out.cobotr_b8 <- coefficients(fit)[9]
out.msath_b9 <- coefficients(fit)[10]
out.fsath_b10 <- coefficients(fit)[11]
out.fhalle_b11 <- coefficients(fit)[12]
out.dmode_b12 <- coefficients(fit)[14]
out.gest_b13 <- coefficients(fit)[15]
out.mage_b14 <- coefficients(fit)[16]
out.SES2_b15<- coefficients(fit)[17]
out.SES3_b16 <- coefficients(fit)[18]
out.SES4_b17 <- coefficients(fit)[19]
out.SES5_b18 <- coefficients(fit)[20]
out.msmok_b19 <- coefficients(fit)[22]



vals <-  c(vals,as.numeric(c(out_b0,out.exp_b1,out.sex_b2,out.sibs1_b3,out.sibs2_b4,out.sibs3_b5,out.LBW_b6,
                             out.cobasi_b7,out.cobotr_b8, out.msath_b9,out.fsath_b10,out.fhalle_b11,
                             out.dmode_b12,out.gest_b13, out.mage_b14,out.SES2_b15,
                             out.SES3_b16,out.SES4_b17,out.SES5_b18,out.msmok_b19)))
para_name <-c(para_name,c("out_b0","out.exp_b1","out.sex_b2","out.sibs1_b3","out.sibs2_b4","out.sibs3_b5","out.LBW_b6",
                          "out.cobasi_b7","out.cobotr_b8", "out.msath_b9","out.fsath_b10","out.fhalle_b11",
                          "out.dmode_b12","out.gest_b13", "out.mage_b14","out.SES2_b15",
                          "out.SES3_b16","out.SES4_b17","out.SES5_b18","out.msmok_b19"))




#response  indicator 

dat$R_y <- ifelse(is.na(dat$asthma_dx_current_6y),0,1)


fit <- glm(R_y ~any_brfeed_an+sexii+number_sibs+LBW_infants+parent_cob_3cat+mother_asthmaii+father_asthmaii+
             family_history_allergyii+
             mode_delivery+pre_term+mother_age_atbirth+as.factor(SEIFA_quint)+
             smoke_mother_pregnancyii,family=binomial(link='logit'),data=dat)


Ry_b0 <- coefficients(fit)[1]
Ry.exp_b1 <- coefficients(fit)[2]
Ry.sex_b2 <- coefficients(fit)[3]
Ry.sibs1_b3 <- coefficients(fit)[4]
Ry.sibs2_b4 <- coefficients(fit)[5]
Ry.sibs3_b5 <- coefficients(fit)[6]
Ry.LBW_b6 <- coefficients(fit)[7]
Ry.cobasi_b7 <- coefficients(fit)[8]
Ry.cobotr_b8 <- coefficients(fit)[9]
Ry.msath_b9 <- coefficients(fit)[10]
Ry.fsath_b10 <- coefficients(fit)[11]
Ry.fhalle_b11 <- coefficients(fit)[12]
Ry.dmode_b12 <- coefficients(fit)[14]
Ry.gest_b13 <- coefficients(fit)[15]
Ry.mage_b14 <- coefficients(fit)[16]
Ry.SES2_b15<- coefficients(fit)[17]
Ry.SES3_b16 <- coefficients(fit)[18]
Ry.SES4_b17 <- coefficients(fit)[19]
Ry.SES5_b18 <- coefficients(fit)[20]
Ry.msmok_b19 <- coefficients(fit)[22]



vals <-  c(vals,as.numeric(c(Ry_b0,Ry.exp_b1,Ry.sex_b2,Ry.sibs1_b3,Ry.sibs2_b4,Ry.sibs3_b5,Ry.LBW_b6,
                             Ry.cobasi_b7,Ry.cobotr_b8, Ry.msath_b9,Ry.fsath_b10,Ry.fhalle_b11,
                             Ry.dmode_b12,Ry.gest_b13, Ry.mage_b14,Ry.SES2_b15,
                             Ry.SES3_b16,Ry.SES4_b17,Ry.SES5_b18,Ry.msmok_b19)))
para_name <-c(para_name,c("Ry_b0","Ry.exp_b1","Ry.sex_b2","Ry.sibs1_b3","Ry.sibs2_b4","Ry.sibs3_b5","Ry.LBW_b6",
                          "Ry.cobasi_b7","Ry.cobotr_b8", "Ry.msath_b9","Ry.fsath_b10","Ry.fhalle_b11",
                          "Ry.dmode_b12","Ry.gest_b13", "Ry.mage_b14","Ry.SES2_b15",
                          "Ry.SES3_b16","Ry.SES4_b17","Ry.SES5_b18","Ry.msmok_b19"))


param <- data.frame(para_name,vals)

param$vals <- round(param$vals,2)

write.csv(param, "parameters.csv")