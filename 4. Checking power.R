##################################################################################
#  Methodology for systematic identification of multiple bias analysis in causal 
#  inference : Checking power and setting effect sizes - 
#  Author: Rushani Wijesuriya 9th of April 2025                           #
##################################################################################

rm(list=ls())


library('fastDummies')
source("gcomp_boot_functions.R")
library(boot)
source("sim.gen.R")


#import parameters 
key_BP <- read.csv("BP_sim.csv")

#-------------------------------realistic scenario

#generate data 

j=2
realistic <- sim.gen(n=2000,nsim=1000,seed=839295,p_U=key_BP[2,j],OR_AU=key_BP[3,j],OR_YU=key_BP[4,j],OR_YAE=key_BP[5,j],p_E=key_BP[6,j],
              sens_X=key_BP[7,j],spec_X=key_BP[8,j], sens_Y=key_BP[9,j],spec_Y=key_BP[10,j],OR_YRy=key_BP[11,j],alpha=key_BP[12,j])


#check power
mylist <- realistic$dataframes_E

#risk difference
mods1 <- lapply(mylist, function(d) {
  glm(glm(asthma_dx_current_6y~bf_true+sexii+as.factor(number_sibs)+mother_age_atbirth+
            LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
            family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
            smoke_mother_pregnancyii+U_gh,family = gaussian("identity"), data=d))
})

p_val_RD <- c()
RD <- c()
for (i in 1:1000){p_val_RD[i] <- summary(mods1[[i]])[["coefficients"]][2,4]
RD[i] <-summary(mods1[[i]])[["coefficients"]][2,1]} 
 

sum(p_val_RD<0.05)/1000
mean(RD)

#risk ratio
mods2 <- lapply(mylist, function(d) {
  glm(asthma_dx_current_6y~bf_true+U_gh+sexii+as.factor(number_sibs)+mother_age_atbirth+
        LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
        family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
        smoke_mother_pregnancyii, data=d, family=binomial(logit))
})



p_val_RR <- c()
RR <- c()
for (i in 1:1000){p_val_RR[i] <- summary(mods2[[i]])[["coefficients"]][2,4]
RR[i] <-summary(mods2[[i]])[["coefficients"]][2,1]} 

sum(p_val_RR<0.05)/1000
mean(exp(RR))

#-------------------------------enhanced scenario

#generate data

j=3
enhanced <- sim.gen(n=2000,nsim=1000,seed=839295,p_U=key_BP[2,j],OR_AU=key_BP[3,j],OR_YU=key_BP[4,j],OR_YAE=key_BP[5,j],p_E=key_BP[6,j],
              sens_X=key_BP[7,j],spec_X=key_BP[8,j], sens_Y=key_BP[9,j],spec_Y=key_BP[10,j],OR_YRy=key_BP[11,j],alpha=key_BP[12,j])


#check power
mylist <- realistic$dataframes_E

#risk difference
mods1 <- lapply(mylist, function(d) {
  glm(glm(asthma_dx_current_6y~bf_true+sexii+as.factor(number_sibs)+mother_age_atbirth+
            LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
            family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
            smoke_mother_pregnancyii+U_gh,family = gaussian("identity"), data=d))
})

p_val_RD <- c()
RD <- c()
for (i in 1:1000){p_val_RD[i] <- summary(mods1[[i]])[["coefficients"]][2,4]
RD[i] <-summary(mods1[[i]])[["coefficients"]][2,1]} 


sum(p_val_RD<0.05)/1000
mean(RD)

#risk ratio
mods2 <- lapply(mylist, function(d) {
  glm(asthma_dx_current_6y~bf_true+U_gh+sexii+as.factor(number_sibs)+mother_age_atbirth+
        LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
        family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
        smoke_mother_pregnancyii, data=d, family=binomial(logit))
})



p_val_RR <- c()
RR <- c()
for (i in 1:1000){p_val_RR[i] <- summary(mods2[[i]])[["coefficients"]][2,4]
RR[i] <-summary(mods2[[i]])[["coefficients"]][2,1]} 

sum(p_val_RR<0.05)/1000
mean(exp(RR))

