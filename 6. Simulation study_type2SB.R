##################################################################################
#  Methodology for systematic identification of multiple bias analysis in causal 
#  inference : simulation study 
#  Data used: HealthNuts data 
#  Author: Rushani Wijesuriya 4th of December 2024                          #
##################################################################################
rm(list=ls())


library('fastDummies')
library(dplyr)
#library(boot)
library(survey)
#library(sandwich)

#set working directory 
#setwd("C:/Users/rushani.wijesuriya/OneDrive - Murdoch Children's Research Institute/QBA-methods paper/Simulation study")

source("sim.gen.R")


#import parameters 
key_BP <- read.csv("BP_sim.csv")

expit <- function(p) {
  exp(p) / (1 + exp(p))
}


RD_SB2 <- list()
RR_SB2 <- list()

scenario <- c("realistic","enhanced")

for (j in 2:3){
  
  #list to save results from each scenario
  RD_BP <- list()
  RR_BP <- list()
  
  for (k in 1:2){
    

#---------------------- generate a large synthetic population to obtain bias parameters----------------

#1. generate synthetic pop (N=1000000)

BP <- sim.gen(n=1000000,nsim=1,seed=839295,p_U=key_BP[2,j],OR_AU=key_BP[3,j],OR_YU=key_BP[4,j],OR_YAE=key_BP[5,j],p_E=key_BP[6,j],
              sens_X=key_BP[7,j],spec_X=key_BP[8,j], sens_Y=key_BP[9,j],spec_Y=key_BP[10,j],OR_YRy=key_BP[11,j],alpha=key_BP[12,j])


##2.Obtain bias parameters

# PS2 model
E_model <- glm(E ~ bf_star + out_star +bf_star* out_star, data = BP$dataframes_E[[1]], family = binomial(link="logit"))


E1_0     <- coef(E_model)[1]*k
E1_xstar <- coef(E_model)[2]*k
E1_ystar <- coef(E_model)[3]*k
E1_xystar<- coef(E_model)[4]*k

#---------------------- Simulations----------------------------------------------------------------------

#1, Generate data 
dat <- sim.gen(n=2000,nsim=1000,seed=9347252,p_U=key_BP[2,j],OR_AU=key_BP[3,j],OR_YU=key_BP[4,j],OR_YAE=key_BP[5,j],p_E=key_BP[6,j],
               sens_X=key_BP[7,j],spec_X=key_BP[8,j], sens_Y=key_BP[9,j],spec_Y=key_BP[10,j],OR_YRy=key_BP[11,j],alpha=key_BP[12,j])


RR_est <-RR_sd <- RR_uci <-RR_lci <-  c()
RD_est <-RD_sd <-RD_uci <- RD_lci <- c()


for (i in 1:1000){
  
  print(i)
  df <- dat$dataframes[[i]]
  n1 <- nrow(df)
  #apply bias-adjustment method
  df <- df %>% 
      mutate( pE = expit(E1_0+E1_xstar*bf_star+E1_ystar*out_star+E1_xystar*(out_star*bf_star))
              
      )
    
    
    df <- df %>% filter(!is.na(df$pE))
    df$ipw <- 1/df$pE
    
    
    #-------------------risk ratio
    model_RR <- svyglm(out_star~bf_star+sexii+as.factor(number_sibs)+mother_age_atbirth+
                         LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                         family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                         smoke_mother_pregnancyii,family = poisson("log"),
                       design=svydesign(id=~1, weights=~ipw, data=df))
    
    RR_est <- c(RR_est,coefficients(model_RR)[2])
    
    #generate SE (no need to use the sandwich estimator as svy already use this)
    RR_sd <-c(RR_sd,summary(model_RR)$coefficients[2,2])
    
    
    
    #generate CI 
    CI <- c((summary(model_RR)$coef[2, 1] + summary(model_RR)$coefficients[2,2] * qnorm(.025)), summary(model_RR)$coef[2, 1] +
              summary(model_RR)$coefficients[2,2] * qnorm(.975))
    
    
    RR_lci <- c(RR_lci,exp(CI)[1])
    RR_uci<- c(RR_uci,exp(CI)[2])
    
    
    
    #risk difference
    model_RD <- svyglm(out_star~bf_star+sexii+as.factor(number_sibs)+mother_age_atbirth+
                         LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                         family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                         smoke_mother_pregnancyii,family = gaussian("identity"),
                       design=svydesign(id=~1, weights=~ipw, data=df))
    
    
    RD_est <- c(RD_est,coefficients(model_RD)[2])
    
    RD_sd <- c(RD_sd,summary(model_RD)$coefficients[2,2])
    
    #generate CI 
    CI <- c((summary(model_RD)$coef[2, 1] + summary(model_RD)$coefficients[2,2]* qnorm(.025)), summary(model_RD)$coef[2, 1] +
              summary(model_RD)$coefficients[2,2] * qnorm(.975))
    
    RD_lci <- c(RD_lci,CI[1])
    RD_uci <- c(RD_uci,CI[2])
    
  } 

RD_results <- data.frame(rbind(RD_est,RD_sd,RD_lci,RD_uci))
RR_results <- data.frame(rbind(RR_est,RR_sd,RR_lci,RR_uci))

colnames(RD_results) <- colnames(RR_results) <- c(1:1000)

RD_results$factor <-RR_results$factor <- k

RD_BP[[k]] <- RD_results
RR_BP[[k]] <- RR_results
  }
  
  
  RD_BP_df <- rbind(RD_BP[[1]],RD_BP[[2]])
  RD_BP_df$scenario <- scenario[j-1]
  RR_BP_df <- rbind(RR_BP[[1]],RR_BP[[2]])
  RR_BP_df$scenario <- scenario[j-1]
  
  
  RD_SB2[[j-1]] <- RD_BP_df 
  RR_SB2[[j-1]] <- RR_BP_df
}


#save results 
final_RD <- rbind(RD_SB2[[1]],RD_SB2[[2]])
final_RR <- rbind(RR_SB2[[1]],RR_SB2[[2]])


write.csv(final_RD,"SB2_RD.csv")
write.csv(final_RR,"SB2_RR.csv")
