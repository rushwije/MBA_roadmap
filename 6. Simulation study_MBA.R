##################################################################################
#  Methodology for systematic identification of multiple bias analysis in causal 
#  inference : simulation study- simultaneous bias adjustment  
#  Data used: HealthNuts data 
#  Author: Rushani Wijesuriya 4th of Decmeber 2024                          #
##################################################################################

rm(list=ls())


library('fastDummies')
library(dplyr)
library(boot)
library(survey)
library(sandwich)

#set working directory 
#setwd("C:/Users/rushani.wijesuriya/OneDrive - Murdoch Children's Research Institute/QBA-methods paper/Simulation study")

expit <- function(p) {
  exp(p) / (1 + exp(p))
}


#import parameters 
key_BP <- read.csv("BP_sim.csv")

source("sim.gen.R")

RD_MB <- list()
RR_MB <- list()

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

E_model <- glm(E ~ bf_star + out_star +bf_star* out_star, data = BP$dataframes_E[[1]], family = binomial(link="logit"))


E1_0     <- coef(E_model)[1]*k
E1_xstar <- coef(E_model)[2]*k
E1_ystar <- coef(E_model)[3]*k
E1_xystar<- coef(E_model)[4]*k

#U model
u_model <- glm(U_gh~bf_true+asthma_dx_current_6y,data = BP$dataframes[[1]], family = binomial(link="logit"))

u1_0 <- coef(u_model)[1]*k
u1_x <- coef(u_model)[2]*k
u1_y <- coef(u_model)[3]*k


#y model 
y_model <- glm(asthma_dx_current_6y~bf_true+out_star+ sexii+as.factor(number_sibs)+mother_age_atbirth+
                 LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                 family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                 smoke_mother_pregnancyii,data = BP$dataframes[[1]], family = binomial(link="logit"))


y_0 <- coef(y_model)[1]*k
y_x <- coef(y_model)[2]*k
y_ystar <- coef(y_model)[3]*k
y_sex <- coef(y_model)[4]*k
y_sibs1 <- coef(y_model)[5]*k
y_sibs2 <- coef(y_model)[6]*k
y_sibs3<- coef(y_model)[7]*k
y_mage<- coef(y_model)[8]*k
y_LBW<- coef(y_model)[9]*k
y_cob1<- coef(y_model)[10]*k
y_cob2<- coef(y_model)[11]*k
y_masthma<- coef(y_model)[12]*k
y_fasthma<- coef(y_model)[13]*k
y_fhxall<- coef(y_model)[14]*k
y_mdel<- coef(y_model)[15]*k
y_gest<- coef(y_model)[16]*k
y_SES1<- coef(y_model)[17]*k
y_SES2<- coef(y_model)[18]*k
y_SES3<- coef(y_model)[19]*k
y_SES4<- coef(y_model)[20]*k
y_msmok<- coef(y_model)[21]*k


x_model <- glm(bf_true~bf_star+out_star+ sexii+as.factor(number_sibs)+mother_age_atbirth+
                 LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                 family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                 smoke_mother_pregnancyii,data = BP$dataframes[[1]], family = binomial(link="logit"))


x_0     <- coef(x_model)[1]*k
x_xstar <- coef(x_model)[2]*k
x_ystar <- coef(x_model)[3]*k
x_sex <- coef(x_model)[4]*k
x_sibs1 <- coef(x_model)[5]*k
x_sibs2 <- coef(x_model)[6]*k
x_sibs3<- coef(x_model)[7]*k
x_mage<- coef(x_model)[8]*k
x_LBW<- coef(x_model)[9]*k
x_cob1<- coef(x_model)[10]*k
x_cob2<- coef(x_model)[11]*k
x_masthma<- coef(x_model)[12]*k
x_fasthma<- coef(x_model)[13]*k
x_fhxall<- coef(x_model)[14]*k
x_mdel<- coef(x_model)[15]*k
x_gest<- coef(x_model)[16]*k
x_SES1<- coef(x_model)[17]*k
x_SES2<- coef(x_model)[18]*k
x_SES3<- coef(x_model)[19]*k
x_SES4<- coef(x_model)[20]*k
x_msmok<- coef(x_model)[21]*k


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

dat <- sim.gen(n=2000,nsim=1000,seed=14082023,p_U=key_BP[2,j],OR_AU=key_BP[3,j],OR_YU=key_BP[4,j],OR_YAE=key_BP[5,j],p_E=key_BP[6,j],
               sens_X=key_BP[7,j],spec_X=key_BP[8,j], sens_Y=key_BP[9,j],spec_Y=key_BP[10,j],OR_YRy=key_BP[11,j],alpha=key_BP[12,j])


RR_est <-RR_sd <- RR_uci <-RR_lci <-  c()
RD_est <-RD_sd <-RD_uci <- RD_lci <- c()


for (i in 1:1000){
  
  print(i)
  df <- dat$dataframes[[i]]
  
  n1 <- nrow(df)
    df <- df %>% 
      mutate( 
        x_pred=rbinom(n1 , 1, expit(x_0+ x_xstar*bf_star+ x_sex*sexii+ x_sibs1*number_sibs_1+x_sibs2*number_sibs_2+
                                      x_sibs3*number_sibs_3+ x_mage*mother_age_atbirth+x_LBW*LBW_infants+
                                      x_cob1*parent_cob_3cat_1+ x_cob2*parent_cob_3cat_2+x_masthma*mother_asthmaii+
                                      x_fasthma*father_asthmaii+ x_fhxall*family_history_allergyii+
                                      x_mdel*mode_delivery+x_gest*pre_term+x_SES1*SEIFA_quint_1+x_SES2*SEIFA_quint_2+
                                      x_SES3*SEIFA_quint_3+x_SES4*SEIFA_quint_4+ x_msmok*smoke_mother_pregnancyii
                                      )),
        y_pred= rbinom(n1 , 1, expit(y_0+ y_x*x_pred+ y_ystar*out_star+y_sex*sexii+y_sibs1*number_sibs_1+y_sibs2*number_sibs_2+
                                       y_sibs3*number_sibs_3+ y_mage*mother_age_atbirth+y_LBW*LBW_infants+
                                       y_cob1*parent_cob_3cat_1+ y_cob2*parent_cob_3cat_2+y_masthma*mother_asthmaii+
                                       y_fasthma*father_asthmaii+ y_fhxall*family_history_allergyii+
                                       y_mdel*mode_delivery+y_gest*pre_term+y_SES1*SEIFA_quint_1+y_SES2*SEIFA_quint_2+
                                       y_SES3*SEIFA_quint_3+y_SES4*SEIFA_quint_4+ y_msmok*smoke_mother_pregnancyii
                                       )),
         u_pred=rbinom(n1 , 1, expit(u1_0+u1_x*x_pred+u1_y*y_pred)),
         pE = expit(E1_0+E1_xstar*bf_star+E1_ystar*out_star+E1_xystar*(out_star*bf_star)),
        pS1 = expit(Ry_0 + Ry_xstar*bf_star+
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
    
    
    
    
    df <- df %>% filter(!is.na(df$pS1))
    df$ipw <- 1/(df$pS1*df$pE)
    
    
    #-------------------risk ratio
    model_RR <- svyglm(y_pred~x_pred+sexii+as.factor(number_sibs)+mother_age_atbirth+
                         LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                         family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                         smoke_mother_pregnancyii+u_pred,family = poisson("log"),
                       design=svydesign(id=~1, weights=~ipw, data=df))
    
    RR_est <- c(RR_est,coefficients(model_RR)[2])
    
    #generate SE (using sandwitch estimator)
    RR_sd <-c(RR_sd,summary(model_RR)$coefficients[2,2])
    
    #generate CI 
    CI <- c((summary(model_RR)$coef[2, 1] + summary(model_RR)$coefficients[2,2] * qnorm(.025)), summary(model_RR)$coef[2, 1] +
              summary(model_RR)$coefficients[2,2] * qnorm(.975))
    
    RR_lci <- c(RR_lci,exp(CI)[1])
    RR_uci<- c(RR_uci,exp(CI)[2])
    
    
    #risk difference
    
    model_RD <- svyglm(y_pred~x_pred+sexii+as.factor(number_sibs)+mother_age_atbirth+
                         LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                         family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                         smoke_mother_pregnancyii+u_pred,family = gaussian("identity"),
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
  
  
  RD_MB[[j-1]] <- RD_BP_df 
  RR_MB[[j-1]] <- RR_BP_df
}


#save results 
final_RD <- rbind(RD_MB[[1]],RD_MB[[2]])
final_RR <- rbind(RR_MB[[1]],RR_MB[[2]])


write.csv(final_RD,"MB_RD.csv")
write.csv(final_RR,"MB_RR.csv")

