setwd("/Users/Neil/Desktop/Clean_Data/")
library(boot)
library(ggpubr)
library(EnvStats)
library(rstatix)
library(psych)
library(reshape2)
library(nlme)
library(lme4)
library(lattice)
library(dplyr)

# load data
data <- read.csv('Data_combined.csv')

# multiply small DAR coefficients by 1000
data$DUO_PLA_DAR <- data$DUO_PLA_DAR * 1000
data$SOLO_PLA_DAR <- data$SOLO_PLA_DAR * 1000
data$DUO_OXY_DAR <- data$DUO_OXY_DAR * 1000
data$SOLO_OXY_DAR <- data$SOLO_OXY_DAR * 1000


# create average CFC variables for OXY and PLA conditions
DAR_OXY_mean <- rowMeans(data[, c('SOLO_OXY_DAR', 'DUO_OXY_DAR')])
DAR_PLA_mean <- rowMeans(data[, c('SOLO_PLA_DAR', 'DUO_PLA_DAR')])
data <- cbind(data, DAR_PLA_mean, DAR_OXY_mean)

AAC_Z_OXY_mean <- rowMeans(data[, c('SOLO_OXY_AAC_Z', 'DUO_OXY_AAC_Z')])
AAC_Z_PLA_mean <- rowMeans(data[, c('SOLO_PLA_AAC_Z', 'DUO_PLA_AAC_Z')])
data <- cbind(data, AAC_Z_OXY_mean, AAC_Z_PLA_mean)

dPAC_Z_OXY_mean <- rowMeans(data[, c('SOLO_OXY_dPAC_Z', 'DUO_OXY_dPAC_Z')])
dPAC_Z_PLA_mean <- rowMeans(data[, c('SOLO_PLA_dPAC_Z', 'DUO_PLA_dPAC_Z')])
data <- cbind(data, dPAC_Z_OXY_mean, dPAC_Z_PLA_mean)


# create overall average CFC score
DAR_mean <- rowMeans(data[, c('DAR_PLA_mean', 'DAR_OXY_mean')])
AAC_Z_mean <- rowMeans(data[, c('AAC_Z_OXY_mean', 'AAC_Z_PLA_mean')])
dPAC_Z_mean <- rowMeans(data[, c('dPAC_Z_OXY_mean', 'dPAC_Z_PLA_mean')])
data <- cbind(data, DAR_mean, AAC_Z_mean, dPAC_Z_mean)

# table of HSA vs LSA individuals
table(data$median_LAS_total)

# descriptive statistics
describe(data[c('age', 'LAS_total', 'DAR_PLA_mean', 'DAR_OXY_mean', 'AAC_Z_PLA_mean', 'AAC_Z_OXY_mean',
                'dPAC_Z_OXY_mean', 'dPAC_Z_PLA_mean')])


###################### CORRELATION ANALYSIS #########################

# check for normality of the variables
shapiro.test(data$LAS_total)
shapiro.test(data$DAR_mean) # not normally distributed
# also in plots
ggqqplot(data$LAS_total, ylab = "LSAS score")
ggqqplot(data$DAR_mean, ylab = "Delta-Beta CFC (DAR method)")

# test for significant correlation
cor.test(data$LAS_total, data$DAR_mean, method = "spearman")

# plot and check for linearity of the relationship
ggscatter(data, x = "LAS_total", y = "DAR_mean", 
          add = "reg.line", conf.int = F, 
          cor.coef = T, cor.method = "spearman",
          xlab = "LSAS score", ylab = "Delta-Beta CFC (DAR method)")

# CORRELATION AAC_Z_mean

# check for normality of the variables
shapiro.test(data$LAS_total)
shapiro.test(data$AAC_Z_mean) # both normally distributed
# also in plots
ggqqplot(data$LAS_total, ylab = "LSAS score")
ggqqplot(data$AAC_Z_mean, ylab = "Delta-Beta CFC (AAC method)")

# test for significant correlation
cor.test(data$LAS_total, data$AAC_Z_mean, method = "spearman")

# plot and check for linearity of the relationship
ggscatter(data, x = "LAS_total", y = "AAC_Z_mean", 
          add = "reg.line", conf.int = F, 
          cor.coef = T, cor.method = "spearman",
          xlab = "LSAS score", ylab = "Delta-Beta CFC (AAC method)")

# CORRELATION dPAC_Z_mean

# check for normality of the variables
shapiro.test(data$LAS_total)
shapiro.test(data$dPAC_Z_mean) # not normally distributed
# also in plots
ggqqplot(data$LAS_total, ylab = "LSAS score")
ggqqplot(data$dPAC_Z_mean, ylab = "Delta-Beta CFC (dPAC method)")

# test for significant correlation
cor.test(data$LAS_total, data$dPAC_Z_mean, method = "spearman")

# plot and check for linearity of the relationship
ggscatter(data, x = "LAS_total", y = "dPAC_Z_mean", 
          add = "reg.line", conf.int = F, 
          cor.coef = T, cor.method = "spearman",
          xlab = "LSAS score", ylab = "Delta-Beta CFC (dPAC method)")




###################### LOGISTIC REGRESSION #########################


# LOGISTIC REGRESSION AAC_Z

# run logistic regression model predicting social anxiety group from AAC_Z
set.seed(1649566)
model_AAC = glm(median_LAS_total ~ AAC_Z_mean, data = data, family = binomial)
summary(model_AAC)

probabilities <- predict(model_AAC, type = "response")
logit_AAC = log(probabilities/(1-probabilities))

data = cbind(data, logit_AAC)


ggplot(data, aes(logit_AAC, AAC_Z_mean)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  theme_bw()
  

plot(model_AAC)

# LOGISTIC REGRESSION dPAC_Z

# run logistic regression model predicting social anxiety group from dPAC_Z
set.seed(1649566)
model_dPAC = glm(median_LAS_total ~ dPAC_Z_mean, data = data, family = binomial)
summary(model_dPAC)

probabilities <- predict(model_dPAC, type = "response")
logit_dPAC = log(probabilities/(1-probabilities))

data <-  cbind(data, logit_dPAC)

ggplot(data, aes(logit_dPAC, dPAC_Z_mean))+
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  theme_bw()

plot(model_dPAC)


# LOGISTIC REGRESSION DAR

# run logistic regression model predicting social anxiety group from DAR
set.seed(1649566)
model_DAR = glm(median_LAS_total ~ DAR_mean, data = data, family = binomial)
summary(model_DAR)

probabilities <- predict(model_DAR, type = "response")
logit_DAR = log(probabilities/(1-probabilities))

data <- cbind(data, logit_DAR)

ggplot(data, aes(logit_DAR, DAR_mean)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  theme_bw()
  

plot(model_DAR)

########################### Multilevel Analysis ###################################

# convert dataset to long format
data$ppn <- factor(data$ppn)
DAR_data_wide <- data[, c("ppn", "age", "order_adm", "LAS_total", "median_LAS_total", "DUO_PLA_DAR", "SOLO_PLA_DAR", "DUO_OXY_DAR", "SOLO_OXY_DAR")]
AAC_data_wide <- data[, c("ppn", "age", "order_adm", "LAS_total", "median_LAS_total", "DUO_PLA_AAC_Z", "SOLO_PLA_AAC_Z", "DUO_OXY_AAC_Z", "SOLO_OXY_AAC_Z")]
dPAC_data_wide <- data[, c("ppn", "age", "order_adm", "LAS_total", "median_LAS_total", "DUO_PLA_dPAC_Z", "SOLO_PLA_dPAC_Z", "DUO_OXY_dPAC_Z", "SOLO_OXY_dPAC_Z")]


# Specify id.vars: the variables to keep but not split apart on
AAC_data_long <- melt(AAC_data_wide, id.vars=c("ppn", "age", "order_adm", "LAS_total", "median_LAS_total"))
dPAC_data_long <- melt(dPAC_data_wide, id.vars=c("ppn", "age", "order_adm", "LAS_total", "median_LAS_total"))
DAR_data_long <- melt(DAR_data_wide, id.vars=c("ppn", "age", "order_adm", "LAS_total", "median_LAS_total"))

data_long <- cbind(AAC_data_long, dPAC_data_long$value, DAR_data_long$value)
colnames(data_long) <- c("ppn", "age", "order_adm", "LAS_total", "median_LAS_total", "condition", "AAC_Z", "dPAC_Z", "DAR")


flanker_cond <- c(rep("DUO", 24), rep("SOLO", 24), rep("DUO", 24), rep("SOLO", 24))
oxy_pla <- c(rep("PLA", 48), rep("OXY", 48))

data_long <- cbind(data_long, oxy_pla, flanker_cond)

str(data_long)


# normality of the outcome variable
hist(data_long$LAS_total, xlab = "LSAS score", main = "Histogram of the LSAS Scores")


# Spaghetti plot per group
g1 <- ggplot(data = data_long, aes(x = flanker_cond, y = AAC_Z, group = ppn)) + geom_line()
g2 <- g1 + stat_summary(fun = mean, aes(group = 1), geom = "line", lwd = 3)
g3 <- g2 + facet_grid(. ~ oxy_pla)
print(g3)

g1 <- ggplot(data = data_long, aes(x = flanker_cond, y = dPAC_Z, group = ppn)) + geom_line()
g2 <- g1 + stat_summary(fun = mean, aes(group = 1), geom = "line", lwd = 3)
g3 <- g2 + facet_grid(. ~ oxy_pla)
print(g3)

g1 <- ggplot(data = data_long, aes(x = flanker_cond, y = DAR, group = ppn)) + geom_line()
g2 <- g1 + stat_summary(fun = mean, aes(group = 1), geom = "line", lwd = 3)
g3 <- g2 + facet_grid(. ~ oxy_pla)
print(g3)


# AAC models
# regular regression with only intercept
AAC_0 = gls(AAC_Z ~ 1, data = data_long, method = "ML", na.action = "na.omit")
# null model (only random intercept per ppn)
AAC_1 = lme( AAC_Z ~ 1, random =  ~1|ppn, method = "ML", data = data_long, na.action = "na.omit" )
# hierarchical model selection
AAC_2 = lme( AAC_Z ~ LAS_total, random = ~1 + LAS_total|ppn, method = "ML", data = data_long, na.action = "na.omit" )
AAC_3 = lme( AAC_Z ~ LAS_total + oxy_pla, random = ~1 + LAS_total|ppn, method = "ML", data = data_long, na.action = "na.omit" )
AAC_4 = lme( AAC_Z ~ LAS_total + oxy_pla + flanker_cond, random = ~1 + LAS_total|ppn, method = "ML", data = data_long, na.action = "na.omit" )
AAC_5 = lme( AAC_Z ~ LAS_total*oxy_pla + LAS_total*flanker_cond, random = ~1 + LAS_total|ppn, method = "ML", data = data_long, na.action = "na.omit" )

# summary of models
summary(AAC_0)
summary(AAC_1)
summary(AAC_2)
summary(AAC_3)
summary(AAC_4)
summary(AAC_5)
anova(AAC_1, AAC_2, AAC_3, AAC_4)

# assumptions AAC Models

# homoscedasticity of the residuals
plot(AAC_4, main = "Homoscedasticity of the Residuals Model 4 (AAC)")

# check for normality of all the residuals
# extract random effects, residuals and fitted values
ran_int <- ranef(AAC_4)[["(Intercept)"]]
ran_slp <- ranef(AAC_4)[["LAS_total"]]
resid = residuals(AAC_4)

# make plots
oldpar1<-par()
par(mfrow=c(1,3))
qqnorm(ran_int, main = "Intercepts Q-Q Plot")
qqline(ran_int)
qqnorm(ran_slp, main = "Slopes Q-Q Plot")
qqline(ran_slp)
qqnorm(resid, main = "Residuals Q-Q Plot")
qqline(resid)


# extract fitted values
fitvalues = fitted(AAC_2)

# make plot of fitted values and residuals
par(mfrow=c(1,1))
plot( fitvalues , resid, main = "Scatterplot of the Residuals", xlab = "Fit Values", ylab = "Residuals" )

# homoscedasticity for level 1-residuals
par(mfrow=c(1,1))
plot(data_long$condition, resid, main = "Homoscedasticity of the Residuals", xlab = "Condition", ylab = "Residuals", type = "p")



# dPAC models
# regular regression with only intercept
dPAC_0 = gls(dPAC_Z ~ 1, data = data_long, method = "ML", na.action = "na.omit")
# null model (only random intercept per ppn)
dPAC_1 = lmer( dPAC_Z ~ 1 + (1|ppn), REML = F, data = data_long, na.action = "na.omit" )
# hierarchical model selection
dPAC_2 = lmer( dPAC_Z ~ LAS_total + (1 + LAS_total|ppn), REML = F, data = data_long, na.action = "na.omit" )
dPAC_3 = lmer( dPAC_Z ~ LAS_total + oxy_pla + (1 + LAS_total|ppn), REML = F, data = data_long, na.action = "na.omit" )
dPAC_4 = lmer( dPAC_Z ~ LAS_total + oxy_pla + flanker_cond + (1 + LAS_total|ppn), REML = F, data = data_long, na.action = "na.omit" )

# summary of models
summary(dPAC_0)
summary(dPAC_1)
summary(dPAC_2)
summary(dPAC_3)
summary(dPAC_4)
anova(dPAC_1, dPAC_2, dPAC_3, dPAC_4)


# homoscedasticity of the residuals
plot(dPAC_4, main = "Homoscedasticity of the Residuals Model 4 (dPAC)", xlab = "Fit Values", ylab = "Residuals")

# check for normality of all the residuals
# extract random effects, residuals and fitted values
ran_int <- ranef(dPAC_4)$ppn[["(Intercept)"]]
ran_slp <- ranef(dPAC_4)$ppn[["LAS_total"]]
resid = residuals(dPAC_4)

# make plots
oldpar1<-par()
par(mfrow=c(1,3))
qqnorm(ran_int, main = "Intercepts Q-Q Plot")
qqline(ran_int)
qqnorm(ran_slp, main = "Slopes Q-Q Plot")
qqline(ran_slp)
qqnorm(resid, main = "Residuals Q-Q Plot")
qqline(resid)


# DAR models
# regular regression with only intercept
DAR_0 = gls(DAR ~ 1, data = data_long, method = "ML", na.action = "na.omit")
# null model (only random intercept per ppn)
DAR_1 = lmer( DAR ~ 1 + (1|ppn), REML = F, data = data_long, na.action = "na.omit" )
# hierarchical model selection
DAR_2 = lmer( DAR ~ LAS_total + (1 + LAS_total|ppn), REML = F, data = data_long, na.action = "na.omit" )
DAR_3 = lmer( DAR ~ LAS_total + oxy_pla + (1 + LAS_total|ppn), REML = F, data = data_long, na.action = "na.omit" )
DAR_4 = lmer( DAR ~ LAS_total + oxy_pla + flanker_cond + (1 + LAS_total|ppn), REML = F, data = data_long, na.action = "na.omit" )

# summary of models
summary(DAR_0)
summary(DAR_1)
summary(DAR_2)
summary(DAR_3)
summary(DAR_4)
anova(DAR_1, DAR_2, DAR_3, DAR_4)


# homoscedasticity of the residuals
plot(DAR_4, main = "Homoscedasticity of the Residuals Model 4 (DAR)", xlab = "Fit Values", ylab = "Residuals")

# check for normality of all the residuals
# extract random effects, residuals and fitted values
ran_int <- ranef(DAR_4)$ppn[["(Intercept)"]]
ran_slp <- ranef(DAR_4)$ppn[["LAS_total"]]
resid = residuals(DAR_4)

# make plots
oldpar1<-par()
par(mfrow=c(1,3))
qqnorm(ran_int, main = "Intercepts Q-Q Plot")
qqline(ran_int)
qqnorm(ran_slp, main = "Slopes Q-Q Plot")
qqline(ran_slp)
qqnorm(resid, main = "Residuals Q-Q Plot")
qqline(resid)

