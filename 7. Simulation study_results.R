##################################################################################
#  Methodology for systematic identification of multiple bias analysis in causal 
#  inference : simulation study-assesing performance and generating tables
#  Data used: Simulated Data
#  Author: Rushani Wijesuriya 20 August 2023                           #
##################################################################################

rm(list=ls())

library(rsimsum)
library(dplyr)
library(ggplot2)
library(writexl)

#set working directory to where results (from scripts numbered 6)are stored 
setwd("~/Simulation study/results")

#import saved results  
temp <- list.files(pattern = glob2rx("*.csv"))

sim.scen <- c("realistic","enhanced")

#RD
temp1 <- temp[grep("RD", temp)]
dup1 <- c("SB1_RD.csv","SB2_RD.csv","conf_RD.csv","measurX_RD.csv","measurY_RD.csv","MB_RD.csv")

temp1 <- temp1[order(match(temp1, dup1))]

#RR
temp2 <- temp[grep("RR", temp)] 
dup2 <- c("SB1_RR.csv","SB2_RR.csv","conf_RR.csv","measurX_RR.csv","measurY_RR.csv","MB_RR.csv")
temp2 <- temp2[order(match(temp2, dup2))]

Method.names <- c("SB1-only","SB2-only","CB-only","MBx-only","MBy-only","All")

true_RD <- c(-0.09,-0.09)
true_RR <- c(log(0.55),log(0.60))


for (i in 1:2) {
  
  
  for(k in 1:2){

RD_est <- c()
RD_sd <- c()

RR_est <- c()
RR_sd <- c()


# Create empty matrix to hold values of the estimates and compute performance
# measures
SIM_RD <- matrix(NA, nrow = length(temp1) * 1000, ncol = 4)
SIM_RR <- matrix(NA, nrow = length(temp1) * 1000, ncol = 4)


for (m in 1:length(temp1)){
  
  print(m)
  
  #----- Risk difference
  data <- read.csv(temp1[m])
  data <- data %>% filter(factor==k,scenario == sim.scen[i])
  data <- data[2:1001]
  RD_est <- c(RD_est, as.numeric(as.vector(as.matrix(data[1,]))))
  RD_sd <- c(RD_sd,as.numeric(as.vector(as.matrix(data[2,]))))

  
  #------Risk ratio
  data2 <- read.csv(temp2[m])
  data2 <- data2 %>% filter(factor==k,scenario == sim.scen[i])
  data2 <- data2[2:1001]
  RR_est <-  c(RR_est,as.numeric(as.vector(as.matrix(data2[1,]))))
  RR_sd <- c(RR_sd,as.numeric(as.vector(as.matrix(data2[2,]))))

}



SIM_RD[, 3] <- RD_est
SIM_RR[, 3] <- RR_est

SIM_RD[, 4] <- RD_sd
SIM_RR[, 4] <- RR_sd

SIM_RD[, 2] <- SIM_RR[, 2] <- rep(Method.names, each = 1000)
SIM_RD[, 1] <- SIM_RR[, 1] <- rep(1:1000, length(temp1))

colnames(SIM_RD) <- colnames(SIM_RR) <- c("dataset", "Method","b", "se")

SIM_RD <- as.data.frame(SIM_RD)
SIM_RR <- as.data.frame(SIM_RR)


SIM_RD$dataset <- as.numeric(SIM_RD$dataset)
SIM_RD$b <- as.numeric(SIM_RD$b)
SIM_RD$se <- as.numeric(SIM_RD$se)

SIM_RR$dataset <- as.numeric(SIM_RR$dataset)
SIM_RR$b <- as.numeric(SIM_RR$b)
SIM_RR$se <- as.numeric(SIM_RR$se)

SIM_RD <- SIM_RD[order(SIM_RD$dataset), ]
SIM_RR <- SIM_RR[order(SIM_RR$dataset), ]

## Compute the performance measures

##------------------- For the RD
S1 <- simsum(data = SIM_RD, estvarname = "b", true = true_RD[i], se = "se", 
             methodvar = "Method", ref = "All",x = TRUE)

#summary(S1)

#------------------For the RR
S2 <- simsum(data = SIM_RR, estvarname = "b", true = true_RR[i], se = "se", 
             methodvar = "Method", ref = "All",x=TRUE)

#summary(S2)


#summarize results 

RD_results <- data.frame(get_data(S1, stats = c("bias"))[4],
                         get_data(S1, stats = c("bias"))[2],
                         ((get_data(S1, stats = c("bias"))[2])/true_RD[i])*100,
                         get_data(S1, stats = c("empse"))[2],
                         get_data(S1, stats = c("modelse"))[2],
                         get_data(S1, stats = c("cover"))[2],
                         get_data(S1, stats = c("becover"))[2]
)




RR_results <- data.frame(get_data(S2, stats = c("bias"))[4],
                         get_data(S2, stats = c("bias"))[2],
                         ((get_data(S2, stats = c("bias"))[2])/true_RR[i])*100,
                         get_data(S2, stats = c("empse"))[2],
                         get_data(S2, stats = c("modelse"))[2],
                         get_data(S2, stats = c("cover"))[2],
                         get_data(S2, stats = c("becover"))[2]
)

colnames(RR_results) <-colnames(RD_results) <-  c("Method","Bias","RBias","Emp SE","Model SE","Coverage","Bias elim. Coverage")

write_xlsx(RD_results,paste0("RD.results.",sim.scen[i],k,".xlsx"))
write_xlsx(RR_results,paste0("RR.results.",sim.scen[i],k,".xlsx"))

}


}