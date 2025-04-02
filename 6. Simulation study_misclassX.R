##################################################################################
#  Methodology for systematic identification of multiple bias analysis in causal 
#  inference : simulation study- one-at-a-time exposure misclassification  
#  Data used: HealthNuts data 
#  Author: Rushani Wijesuriya 3rd of December 2024                        #
##################################################################################


library('fastDummies')
library(dplyr)
library(boot)
library(survey)
library(sandwich)

#set working directory 
setwd("C:/Users/rushani.wijesuriya/OneDrive - Murdoch Children's Research Institute/QBA-methods paper/Simulation study")

source("sim.gen.R")

#import parameters 
key_BP <- read.csv("BP_sim.csv")


expit <- function(p) {
  exp(p) / (1 + exp(p))
}


RD_measurX <- list()
RR_measurX <- list()

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


#x model
x_model <- glm(bf_true~bf_star+ sexii+as.factor(number_sibs)+mother_age_atbirth+
                 LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                 family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                 smoke_mother_pregnancyii+out_star,data = BP$dataframes[[1]], family = binomial(link="logit"))


x_0     <- coef(x_model)[1]*k
x_xstar <- coef(x_model)[2]*k
x_sex <- coef(x_model)[3]*k
x_sibs1 <- coef(x_model)[4]*k
x_sibs2 <- coef(x_model)[5]*k
x_sibs3<- coef(x_model)[6]*k
x_mage<- coef(x_model)[7]*k
x_LBW<- coef(x_model)[8]*k
x_cob1<- coef(x_model)[9]*k
x_cob2<- coef(x_model)[10]*k
x_masthma<- coef(x_model)[11]*k
x_fasthma<- coef(x_model)[12]*k
x_fhxall<- coef(x_model)[13]*k
x_mdel<- coef(x_model)[14]*k
x_gest<- coef(x_model)[15]*k
x_SES1<- coef(x_model)[16]*k
x_SES2<- coef(x_model)[17]*k
x_SES3<- coef(x_model)[18]*k
x_SES4<- coef(x_model)[19]*k
x_msmok<- coef(x_model)[20]*k
x_y <- coef(x_model)[21]*k


#---------------------- Simulations----------------------------------------------------------------------

#1, Generate data 

dat <- sim.gen(n=2000,nsim=1000,seed=6533800,p_U=key_BP[2,j],OR_AU=key_BP[3,j],OR_YU=key_BP[4,j],OR_YAE=key_BP[5,j],p_E=key_BP[6,j],
               sens_X=key_BP[7,j],spec_X=key_BP[8,j], sens_Y=key_BP[9,j],spec_Y=key_BP[10,j],OR_YRy=key_BP[11,j],alpha=key_BP[12,j])


RR_est <-RR_sd <- RR_uci <-RR_lci <-  c()
RD_est <-RD_sd <-RD_uci <- RD_lci <- c()


#2. Apply bias-adjustment method

for (i in 1:1000){
  
  print(i)
  
  n1 <- nrow(dat$dataframes[[i]])
  
  df <- dat$dataframes[[i]] 
    df <- df %>% 
      mutate( 
        x_pred=rbinom(n1 , 1, expit(x_0+ x_xstar*bf_star+ x_sex*sexii+ x_sibs1*number_sibs_1+x_sibs2*number_sibs_2+
                                       x_sibs3*number_sibs_3+ x_mage*mother_age_atbirth+x_LBW*LBW_infants+
                                       x_cob1*parent_cob_3cat_1+ x_cob2*parent_cob_3cat_2+x_masthma*mother_asthmaii+
                                       x_fasthma*father_asthmaii+ x_fhxall*family_history_allergyii+
                                       x_mdel*mode_delivery+x_gest*pre_term+x_SES1*SEIFA_quint_1+x_SES2*SEIFA_quint_2+
                                       x_SES3*SEIFA_quint_3+x_SES4*SEIFA_quint_4+ x_msmok*smoke_mother_pregnancyii+
                                       x_y*out_star))
        
      )
       
    
    

    #-------------------Fit the analysis model : risk ratio
    model_RR <- glm(out_star~x_pred+sexii+as.factor(number_sibs)+mother_age_atbirth+
                      LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                      family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                      smoke_mother_pregnancyii,family = poisson("log"),data=df)
    
    RR_est <- c(RR_est,coefficients(model_RR)[2])
    
    #generate SE (using sandwitch estimator)
    RR_sd <-c(RR_sd,sqrt(sandwich(model_RR)[2,2]))
    
    
    #generate CI 
    CI <- c((summary(model_RR)$coef[2, 1] + sqrt(sandwich(model_RR)[2,2]) * qnorm(.025)), summary(model_RR)$coef[2, 1] +
              sqrt(sandwich(model_RR)[2,2]) * qnorm(.975))
    
    
    RR_lci <- c(RR_lci,exp(CI)[1])
    RR_uci<- c(RR_uci,exp(CI)[2])
    
    
    #-------------------Fit the analysis model : risk difference 
    
    model_RD <- glm(out_star~x_pred+sexii+as.factor(number_sibs)+mother_age_atbirth+
                      LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                      family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                      smoke_mother_pregnancyii,family = gaussian("identity"), data=df)
    
    RD_est <- c(RD_est,coefficients(model_RD)[2])
    
    RD_sd <- c(RD_sd,sqrt(sandwich(model_RD)[2,2]))
    
    #generate CI 
    CI <- c((summary(model_RD)$coef[2, 1] + sqrt(sandwich(model_RD)[2,2]) * qnorm(.025)), summary(model_RD)$coef[2, 1] +
              sqrt(sandwich(model_RD)[2,2]) * qnorm(.975))
    
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
  
  
  RD_measurX[[j-1]] <- RD_BP_df 
  RR_measurX[[j-1]] <- RR_BP_df
}


#save results 
final_RD <- rbind(RD_measurX[[1]],RD_measurX[[2]])
final_RR <- rbind(RR_measurX[[1]],RR_measurX[[2]])

write.csv(final_RD,"measurX_RD.csv")
write.csv(final_RR,"measurX_RR.csv")


