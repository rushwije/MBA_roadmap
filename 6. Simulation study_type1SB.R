##################################################################################
#  Methodology for systematic identification of multiple bias analysis in causal 
#  inference : simulation study 
#  Data used: HealthNuts data 
#  Author: Rushani Wijesuriya 4th of December 2024                     #
##################################################################################



library('fastDummies')
library(dplyr)
library(boot)
library(survey)
library(sandwich)

#set working directory 
setwd("C:/Users/rushani.wijesuriya/OneDrive - Murdoch Children's Research Institute/QBA-methods paper/Simulation study")

source("sim.gen.R")

expit <- function(p) {
  exp(p) / (1 + exp(p))
}


RD_SB1 <- list()
RR_SB1 <- list()

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

#Ry model 
Ry_model <- glm(Ry ~ bf_star+sexii+as.factor(number_sibs)+mother_age_atbirth+
                 LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                 family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                 smoke_mother_pregnancyii+out_star_comp, data = BP$dataframes[[1]], family = binomial(link="logit"))


Ry_0 <- coef(Ry_model)["(Intercept)"]*k
Ry_xstar <- coef(Ry_model)["bf_star"]*k
Ry_sex <- coef(Ry_model)["sexii"]*k
Ry_sibs1 <- coef(Ry_model)["as.factor(number_sibs)1"]*k
Ry_sibs2 <- coef(Ry_model)["as.factor(number_sibs)2"]*k
Ry_sibs3<- coef(Ry_model)["as.factor(number_sibs)3"]*k
Ry_mage<- coef(Ry_model)["mother_age_atbirth"]*k
Ry_LBW<- coef(Ry_model)["LBW_infants"]*k
Ry_cob1<- coef(Ry_model)["as.factor(parent_cob_3cat)1"]*k
Ry_cob2<- coef(Ry_model)["as.factor(parent_cob_3cat)2"]*k
Ry_masthma<- coef(Ry_model)["mother_asthmaii"]*k
Ry_fasthma<- coef(Ry_model)["father_asthmaii"]*k
Ry_fhxall<- coef(Ry_model)["family_history_allergyii"]*k
Ry_mdel<- coef(Ry_model)["mode_delivery"]*k
Ry_gest<- coef(Ry_model)["pre_term"]*k
Ry_SES1<- coef(Ry_model)["as.factor(SEIFA_quint)1"]*k
Ry_SES2<- coef(Ry_model)["as.factor(SEIFA_quint)2"]*k
Ry_SES3<- coef(Ry_model)["as.factor(SEIFA_quint)3"]*k
Ry_SES4<- coef(Ry_model)["as.factor(SEIFA_quint)4"]*k
Ry_msmok<- coef(Ry_model)["smoke_mother_pregnancyii"]*k
Ry_ystar <- coef(Ry_model)["out_star_comp"]*k

#---------------------- Simulations----------------------------------------------------------------------

#1, Generate data 

dat <- sim.gen(n=2000,nsim=1000,seed=13082023,p_U=key_BP[2,j],OR_AU=key_BP[3,j],OR_YU=key_BP[4,j],OR_YAE=key_BP[5,j],p_E=key_BP[6,j],
               sens_X=key_BP[7,j],spec_X=key_BP[8,j], sens_Y=key_BP[9,j],spec_Y=key_BP[10,j],OR_YRy=key_BP[11,j],alpha=key_BP[12,j])

RR_est <-RR_sd <- RR_uci <-RR_lci <-  c()
RD_est <-RD_sd <-RD_uci <- RD_lci <- c()


for (i in 1:1000){
    print(i)
    df <- dat$dataframes[[i]] 
    n1 <- nrow(df)
   
    #apply bias-adjustment method
    df <-  df %>% mutate( pS = expit(Ry_xstar*bf_star+
                                      Ry_sex*sexii+Ry_sibs1*number_sibs_1+
                                       Ry_sibs2*number_sibs_2+
                                        Ry_sibs3*number_sibs_3+
                                           Ry_mage*mother_age_atbirth+
                                                      Ry_LBW*LBW_infants+
                                                      Ry_cob1*parent_cob_3cat_1+
                                                      Ry_cob2*parent_cob_3cat_2+
                                                      Ry_masthma*mother_asthmaii+
                                                      Ry_fasthma*father_asthmaii+
                                                      Ry_fhxall*family_history_allergyii+
                                                      Ry_mdel*mode_delivery+
                                                      Ry_gest*pre_term+
                                                      Ry_SES1*SEIFA_quint_1+
                                                      Ry_SES2*SEIFA_quint_2+
                                                      Ry_SES3*SEIFA_quint_3+
                                                      Ry_SES4*SEIFA_quint_4+
                                                      Ry_msmok*smoke_mother_pregnancyii+
                                       Ry_ystar*out_star)
                          
    )
    
    
    df <- df %>% filter(!is.na(df$pS))
    df$ipw <- 1/df$pS
    
    
    #-------------------Fit the analysis model : risk ratio
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
    
    
    #-------------------Fit the analysis model : risk difference 
    
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
  
  
  RD_SB1[[j-1]] <- RD_BP_df 
  RR_SB1[[j-1]] <- RR_BP_df
}

#save results 
final_RD <- rbind(RD_SB1[[1]],RD_SB1[[2]])
final_RR <- rbind(RR_SB1[[1]],RR_SB1[[2]])


write.csv(final_RD,"SB1_RD.csv")
write.csv(final_RR,"SB1_RR.csv")



