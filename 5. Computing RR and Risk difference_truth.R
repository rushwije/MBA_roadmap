##################################################################################
#  Methodology for systematic identification of multiple bias analysis in causal 
#  inference : Estimating  the "true" risk difference and the risk ratio
#  Data used: HealthNuts data 
#  Author: Rushani Wijesuriya 17th of August 2023                          #
##################################################################################
library(boot)


expit <- function(p) {
  exp(p) / (1 + exp(p))
}


#import parameters 
key_BP <- read.csv("BP_sim.csv")

source("sim.gen.R")

#step 1: generate large synthetic pop (set j=2 for realistic and j=3 for enhanced)
j=2

synthetic <- sim.gen(n=1000000,nsim=1,seed=13082023,p_U=key_BP[2,j],OR_AU=key_BP[3,j],OR_YU=key_BP[4,j],OR_YAE=key_BP[5,j],p_E=key_BP[6,j],
              sens_X=key_BP[7,j],spec_X=key_BP[8,j], sens_Y=key_BP[9,j],spec_Y=key_BP[10,j],OR_YRy=key_BP[11,j],alpha=key_BP[12,j])



final_diff <- glm(asthma_dx_current_6y~bf_true+U_gh+sexii+as.factor(number_sibs)+mother_age_atbirth+
                       LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                       family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                       smoke_mother_pregnancyii,family=gaussian("identity"), data=synthetic$dataframes_E[[1]])


summary(final_diff)$coeff


final_RR <- glm(asthma_dx_current_6y~bf_true+U_gh+sexii+as.factor(number_sibs)+mother_age_atbirth+
                    LBW_infants+as.factor(parent_cob_3cat)+ mother_asthmaii+father_asthmaii+
                    family_history_allergyii+mode_delivery+pre_term+as.factor(SEIFA_quint)+
                    smoke_mother_pregnancyii,family=poisson("log"), data=synthetic$dataframes_E[[1]])


exp(coef(final_RR))


#realistic
#RD= -0.09
#RR=0.55



#enhanced
#RD= -0.09
#RR=0.60



