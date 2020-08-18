setwd("/Users/Neil/Google Drive/Uni Leiden/MA/Thesis/Analysis/Simulation_study")

library(boot)
library(EnvStats)
library(psych)
library(knitr)


set.seed(1649566)

# import simulated data with PAC
PAC_AAC_P <- read.table('/Users/Neil/Google Drive/Uni Leiden/MA/Thesis/Analysis/Simulation_Study/Results_PAC/PAC_AAC_P_128Hz.txt', sep = ',', header = F,)[, 1]
PAC_AAC_R <- read.table('/Users/Neil/Google Drive/Uni Leiden/MA/Thesis/Analysis/Simulation_Study/Results_PAC/PAC_AAC_R_128Hz.txt', sep = ',', header = F,)[, 1]
PAC_AAC_Z <- read.table('/Users/Neil/Google Drive/Uni Leiden/MA/Thesis/Analysis/Simulation_Study/Results_PAC/PAC_AAC_Z_128Hz.txt', sep = ',', header = F,)[, 1]
PAC_dPAC <- read.table('/Users/Neil/Google Drive/Uni Leiden/MA/Thesis/Analysis/Simulation_Study/Results_PAC/PAC_dPAC_128Hz.txt', sep = ',', header = F,)[, 1]
PAC_dPAC_Z <- read.table('/Users/Neil/Google Drive/Uni Leiden/MA/Thesis/Analysis/Simulation_Study/Results_PAC/PAC_z_128Hz.txt', sep = ',', header = F,)[, 1]

subj <- 1:20
dat_pac <- cbind(subj, PAC_AAC_P, PAC_AAC_R, PAC_AAC_Z, PAC_dPAC, PAC_dPAC_Z)


# import simulated data without PAC
No_PAC_AAC_P <- read.table('/Users/Neil/Google Drive/Uni Leiden/MA/Thesis/Analysis/Simulation_Study/Results_noPAC/NoPAC_AAC_P_128Hz.txt', sep = ',', header = F,)[, 1]
No_PAC_AAC_R <- read.table('/Users/Neil/Google Drive/Uni Leiden/MA/Thesis/Analysis/Simulation_Study/Results_noPAC/NoPAC_AAC_R_128Hz.txt', sep = ',', header = F,)[, 1]
No_PAC_AAC_Z <- read.table('/Users/Neil/Google Drive/Uni Leiden/MA/Thesis/Analysis/Simulation_Study/Results_noPAC/NoPAC_AAC_Z_128Hz.txt', sep = ',', header = F,)[, 1]
No_PAC_dPAC <- read.table('/Users/Neil/Google Drive/Uni Leiden/MA/Thesis/Analysis/Simulation_Study/Results_noPAC/NoPAC_dPAC_128Hz.txt', sep = ',', header = F,)[, 1]
No_PAC_dPAC_Z <- read.table('/Users/Neil/Google Drive/Uni Leiden/MA/Thesis/Analysis/Simulation_Study/Results_noPAC/NoPAC_z_128Hz.txt', sep = ',', header = F,)[, 1]

data <- cbind(dat_pac, No_PAC_AAC_P, No_PAC_AAC_R, No_PAC_AAC_Z, No_PAC_dPAC, No_PAC_dPAC_Z)



# import DAR results
DAR_PAC1 <- read.csv("Results_DAR_PAC/SimulatedData1_PAC_DAR.csv", header = F)
DAR_PAC2 <- read.csv("Results_DAR_PAC/SimulatedData2_PAC_DAR.csv", header = F)
DAR_PAC3 <- read.csv("Results_DAR_PAC/SimulatedData3_PAC_DAR.csv", header = F)
DAR_PAC4 <- read.csv("Results_DAR_PAC/SimulatedData4_PAC_DAR.csv", header = F)
DAR_PAC5 <- read.csv("Results_DAR_PAC/SimulatedData5_PAC_DAR.csv", header = F)
DAR_PAC6 <- read.csv("Results_DAR_PAC/SimulatedData6_PAC_DAR.csv", header = F)
DAR_PAC7 <- read.csv("Results_DAR_PAC/SimulatedData7_PAC_DAR.csv", header = F)
DAR_PAC8 <- read.csv("Results_DAR_PAC/SimulatedData8_PAC_DAR.csv", header = F)
DAR_PAC9 <- read.csv("Results_DAR_PAC/SimulatedData9_PAC_DAR.csv", header = F)
DAR_PAC10 <- read.csv("Results_DAR_PAC/SimulatedData10_PAC_DAR.csv", header = F)
DAR_PAC11 <- read.csv("Results_DAR_PAC/SimulatedData11_PAC_DAR.csv", header = F)
DAR_PAC12 <- read.csv("Results_DAR_PAC/SimulatedData12_PAC_DAR.csv", header = F)
DAR_PAC13 <- read.csv("Results_DAR_PAC/SimulatedData13_PAC_DAR.csv", header = F)
DAR_PAC14 <- read.csv("Results_DAR_PAC/SimulatedData14_PAC_DAR.csv", header = F)
DAR_PAC15 <- read.csv("Results_DAR_PAC/SimulatedData15_PAC_DAR.csv", header = F)
DAR_PAC16 <- read.csv("Results_DAR_PAC/SimulatedData16_PAC_DAR.csv", header = F)
DAR_PAC17 <- read.csv("Results_DAR_PAC/SimulatedData17_PAC_DAR.csv", header = F)
DAR_PAC18 <- read.csv("Results_DAR_PAC/SimulatedData18_PAC_DAR.csv", header = F)
DAR_PAC19 <- read.csv("Results_DAR_PAC/SimulatedData19_PAC_DAR.csv", header = F)
DAR_PAC20 <- read.csv("Results_DAR_PAC/SimulatedData20_PAC_DAR.csv", header = F)

dar_pac1 <- mean(as.matrix(DAR_PAC1))
dar_pac2 <- mean(as.matrix(DAR_PAC2))
dar_pac3 <- mean(as.matrix(DAR_PAC3))
dar_pac4 <- mean(as.matrix(DAR_PAC4))
dar_pac5 <- mean(as.matrix(DAR_PAC5))
dar_pac6 <- mean(as.matrix(DAR_PAC6))
dar_pac7 <- mean(as.matrix(DAR_PAC7))
dar_pac8 <- mean(as.matrix(DAR_PAC8))
dar_pac9 <- mean(as.matrix(DAR_PAC9))
dar_pac10 <- mean(as.matrix(DAR_PAC10))
dar_pac11 <- mean(as.matrix(DAR_PAC11))
dar_pac12 <- mean(as.matrix(DAR_PAC12))
dar_pac13 <- mean(as.matrix(DAR_PAC13))
dar_pac14 <- mean(as.matrix(DAR_PAC14))
dar_pac15 <- mean(as.matrix(DAR_PAC15))
dar_pac16 <- mean(as.matrix(DAR_PAC16))
dar_pac17 <- mean(as.matrix(DAR_PAC17))
dar_pac18 <- mean(as.matrix(DAR_PAC18))
dar_pac19 <- mean(as.matrix(DAR_PAC19))
dar_pac20 <- mean(as.matrix(DAR_PAC20))

DAR1 <- read.csv("Results_DAR_noPAC/SimulatedData1_DAR.csv", header = F)
DAR2 <- read.csv("Results_DAR_noPAC/SimulatedData2_DAR.csv", header = F)
DAR3 <- read.csv("Results_DAR_noPAC/SimulatedData3_DAR.csv", header = F)
DAR4 <- read.csv("Results_DAR_noPAC/SimulatedData4_DAR.csv", header = F)
DAR5 <- read.csv("Results_DAR_noPAC/SimulatedData5_DAR.csv", header = F)
DAR6 <- read.csv("Results_DAR_noPAC/SimulatedData6_DAR.csv", header = F)
DAR7 <- read.csv("Results_DAR_noPAC/SimulatedData7_DAR.csv", header = F)
DAR8 <- read.csv("Results_DAR_noPAC/SimulatedData8_DAR.csv", header = F)
DAR9 <- read.csv("Results_DAR_noPAC/SimulatedData9_DAR.csv", header = F)
DAR10 <- read.csv("Results_DAR_noPAC/SimulatedData10_DAR.csv", header = F)
DAR11 <- read.csv("Results_DAR_noPAC/SimulatedData11_DAR.csv", header = F)
DAR12 <- read.csv("Results_DAR_noPAC/SimulatedData12_DAR.csv", header = F)
DAR13 <- read.csv("Results_DAR_noPAC/SimulatedData13_DAR.csv", header = F)
DAR14 <- read.csv("Results_DAR_noPAC/SimulatedData14_DAR.csv", header = F)
DAR15 <- read.csv("Results_DAR_noPAC/SimulatedData15_DAR.csv", header = F)
DAR16 <- read.csv("Results_DAR_noPAC/SimulatedData16_DAR.csv", header = F)
DAR17 <- read.csv("Results_DAR_noPAC/SimulatedData17_DAR.csv", header = F)
DAR18 <- read.csv("Results_DAR_noPAC/SimulatedData18_DAR.csv", header = F)
DAR19 <- read.csv("Results_DAR_noPAC/SimulatedData19_DAR.csv", header = F)
DAR20 <- read.csv("Results_DAR_noPAC/SimulatedData20_DAR.csv", header = F)

dar1 <- mean(as.matrix(DAR1))
dar2 <- mean(as.matrix(DAR2))
dar3 <- mean(as.matrix(DAR3))
dar4 <- mean(as.matrix(DAR4))
dar5 <- mean(as.matrix(DAR5))
dar6 <- mean(as.matrix(DAR6))
dar7 <- mean(as.matrix(DAR7))
dar8 <- mean(as.matrix(DAR8))
dar9 <- mean(as.matrix(DAR9))
dar10 <- mean(as.matrix(DAR10))
dar11 <- mean(as.matrix(DAR11))
dar12 <- mean(as.matrix(DAR12))
dar13 <- mean(as.matrix(DAR13))
dar14 <- mean(as.matrix(DAR14))
dar15 <- mean(as.matrix(DAR15))
dar16 <- mean(as.matrix(DAR16))
dar17 <- mean(as.matrix(DAR17))
dar18 <- mean(as.matrix(DAR18))
dar19 <- mean(as.matrix(DAR19))
dar20 <- mean(as.matrix(DAR20))

DAR_PAC_mean <- rbind(dar_pac1, dar_pac2, dar_pac3, dar_pac4, dar_pac5,dar_pac6, dar_pac7, dar_pac8,
                             dar_pac9, dar_pac10, dar_pac11, dar_pac12, dar_pac13, dar_pac14, dar_pac15,
                             dar_pac16, dar_pac17, dar_pac18, dar_pac19, dar_pac20)

DAR_No_PAC_mean <- rbind(dar1, dar2, dar3, dar4, dar5, dar6, dar7, dar8, dar9, dar10, dar11,
                                dar12, dar13, dar14, dar15, dar16, dar17, dar18, dar19, dar20)

data <- as.data.frame(cbind(data, DAR_PAC_mean, DAR_No_PAC_mean))
colnames(data)[12:13] <- c('DAR_PAC_mean', 'DAR_No_PAC_mean')


describe_data <- describe(data[c("PAC_AAC_Z", "No_PAC_AAC_Z", "PAC_dPAC_Z", "No_PAC_dPAC_Z", "DAR_PAC_mean", "DAR_No_PAC_mean")])

describe_data[,c(3:5, 8:13)]

# prepare data for logistic regression
condition <- c(rep.int(1, 20), rep.int(0, 20))
DAR <- c(data$DAR_PAC_mean, data$DAR_No_PAC_mean)
dPAC <- c(data$PAC_dPAC_Z, data$No_PAC_dPAC_Z)
AAC <- c(data$PAC_AAC_Z, data$No_PAC_AAC_Z)
CFC_cond <- as.data.frame(cbind(condition, DAR, dPAC, AAC))

# predict condition based on AAC_Z and dPAC_Z
set.seed(1649566)

model_AAC <-  glm(condition ~ AAC, data = CFC_cond)
summary(model_AAC)
model_AAC_cv = cv.glm(CFC_cond, model_AAC , K = 10)
# misclassification rate in training set
model_AAC_cv$delta[1]

model_dPAC <-  glm(condition ~ dPAC, data = CFC_cond)
summary(model_dPAC)
model_dPAC_cv = cv.glm(CFC_cond, model_dPAC , K = 10)
# misclassification rate in training set
model_dPAC_cv$delta[1]

model_DAR = glm(condition ~ DAR, data = CFC_cond, family = binomial)
summary(model_DAR)
model_DAR_cv = cv.glm(CFC_cond, model_DAR , K = 10)
model_DAR_cv
# misclassification rate in training set
model_DAR_cv$delta[1]



# split data in to both conditions
PAC <- sim_dat[which(sim_dat$condition == 1), ]
noPAC <- sim_dat[which(sim_dat$condition == 0), ]

# perform permutation tests
permu_PAC <- twoSamplePermutationTestLocation(data$PAC_dPAC_Z, data$No_PAC_dPAC_Z, n.permutations = 10000, fcn="mean")
permu_PAC

permu_AAC <- twoSamplePermutationTestLocation(data$PAC_AAC_Z, data$No_PAC_AAC_Z, n.permutations = 10000, fcn="mean")
permu_AAC

permu_DAR <- twoSamplePermutationTestLocation(data$DAR_PAC_mean, data$DAR_No_PAC_mean, n.permutations = 10000, fcn="mean")
permu_DAR


boxplot(data[c('PAC_dPAC_Z', 'No_PAC_dPAC_Z')])
boxplot(data[c('PAC_AAC_Z', 'No_PAC_AAC_Z')])
boxplot(data[c('DAR_PAC_mean', 'DAR_No_PAC_mean')])



