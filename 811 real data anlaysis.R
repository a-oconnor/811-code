### Simulation for inverse probability weighting
#Packages=============================================
rm(list=ls())
library(locfit)
library(lme4)
library(dplyr)
library(xtable)
library(survival)
source("functions.R")
#Import code===============================================================================================
#import data
realdat <- read.csv("cleaned_dat.csv")
realdat <- realdat[realdat$feedlot!=5,] # all feedlots satisfies the positivity cond
realdat <- realdat[order(realdat$feedlot),]



#Outcome: Multivariate model===============================================================================================
M1.1 <- glmer(brd50 ~ BVDV_binary + age + weight + mix_final_v1 + BVDV_PI_origin_FINAL_v1 + PV_pre2_final + (1|feedlot),
              family = binomial, data=realdat, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e6)))


# Treatment modelling & weight calculation
#Exposure: Random effects===============================================================================================
# IPW -- clusters are treated as random effects
M2.1 <- glmer(BVDV_binary ~ age + weight + mix_final_v1 + BVDV_PI_origin_FINAL_v1 + PV_pre2_final + (1|feedlot),
              family = binomial, data=realdat, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e6)))
ps.re <- predict(M2.1, type="response")

#Exposure: Fixed effects===============================================================================================
# IPW weighting logistic regression (clusters are treated as fixed effects); weights are the same as the results from glmer
M2.2 <- glm(BVDV_binary ~ age + weight + mix_final_v1 + BVDV_PI_origin_FINAL_v1 + PV_pre2_final + strata(feedlot), 
            family = binomial, data=realdat)
ps.fe <- M2.2$fitted.values

#Combined the dataset==========================
realdat_new <- data.frame(realdat$AnimalID, realdat$feedlot, realdat$cohort, realdat$brd50, realdat$BVDV_binary, realdat$age,
                          realdat$weight, realdat$mix_final_v1, realdat$BVDV_PI_origin_FINAL_v1, realdat$PV_pre2_final, ps.re, ps.fe)
colnames(realdat_new) <- c("AnimalID", "feedlot", "cohort", "Y", "A",
                           "age", "weight", "mix_hist","PIgroup", "PV_vacc", "ps.re",
                           "ps.fe")

write.csv(realdat_new, "feedlotONLY_data_with_weight.csv")

# Outcome model prediction ============================================
Xmat2<- model.matrix(M1.1, data=realdat_new)
head(Xmat2)
tail(Xmat2)
Xmat2[,2]<-1
M1Xmat2<- Xmat2
Xmat2[,2]<-0
M0Xmat2<- Xmat2

M1.2 <- glm(Y ~  0 + strata(feedlot), family = binomial, data=realdat_new)
M1_M1.1 <- expit(M1Xmat2%*% fixef(M1.1) + as.vector(model.matrix(M1.2)%*% as.matrix(ranef(M1.1)$`feedlot`))) # BLUP
M0_M1.1 <- expit(M0Xmat2%*% fixef(M1.1) + as.vector(model.matrix(M1.2)%*% as.matrix(ranef(M1.1)$`feedlot`)))

realdat_new <- data.frame(realdat_new, M0_M1.1, M1_M1.1)

#Risks from Exposure IPW  model===============================================================================================
ipw.re.p0 <- with(realdat_new, mean(Y*(1-A)/(1-ps.re)))
ipw.re.p1 <- with(realdat_new, mean(Y*A/ps.re))
ipw.fe.p0 <- with(realdat_new, mean(Y*(1-A)/(1-ps.fe)))
ipw.fe.p1 <- with(realdat_new, mean(Y*A/ps.fe))

#Risks from Direct/Univariate Outcome model===============================================================================================
# model based est
realdat_new_sub1 <- realdat_new[realdat_new$A==1,]
realdat_new_sub0 <- realdat_new[realdat_new$A==0,]
direct.p0 <- with(realdat_new_sub0, mean(Y*(1-A)))
direct.p1 <- with(realdat_new_sub1, mean(Y*A))

#Risks from Model based/ Multivariate  model===============================================================================================

# mean of outcome Model based prediction
M1.1_phat1 <- with(realdat_new_sub1, mean(M1_M1.1))
M1.1_phat0 <- with(realdat_new_sub0, mean(M0_M1.1))

#Create a table ==============================================================================================
p0_vec <- rbind(ipw.re.p0, ipw.fe.p0, 
                direct.p0, M1.1_phat0) 

p1_vec <- rbind(ipw.re.p1, ipw.fe.p1, 
                direct.p1, M1.1_phat1)

tau_vec <- p1_vec - p0_vec
OR_vec <- cbind(OR_func(p0 = p0_vec, p1 = p1_vec))
output1 <- round(cbind(p0_vec, p1_vec, tau_vec, OR_vec),digits=6)

rownames(output1) <- c("IPW random", "IPW fixed", "direct/Univariate", "Mbased/Multivariate"
                      )
output1 <- data.frame(output1)
colnames(output1) <- c("p0", "p1", "p1-p0", "OR")

sink("causal_effect_estimate.txt", append = TRUE)
cat("=================================================\n")
cat("having only feedlot as random/fixed effect \n")
cat("size=", nrow(realdat), "\n")
print(output1)
sink()

# sink("feedlot_only_model_parameter_estimate.txt", append = TRUE)
# cat("=================================================\n")
# cat("Outcome model \n")
# print(summary(M1.1))
# cat("=================================================\n")
# cat("Treatment model random feedlot \n")
# print(summary(M2.1))
# cat("=================================================\n")
# cat("Treatment model fixed feedlot \n")
# print(summary(M2.2))
# cat("=================================================\n")
# cat("conditional model\n")
# print(summary(M2.3))
# sink()
# 
# sink("causal_effect_estimate.txt", append = TRUE)
# cat("=================================================\n")
# cat("Only having feedlot as random/fixed effect \n")
# print(output1)
# sink()

subdat_0 <- realdat_new[realdat_new$A==0,]
subdat_1 <- realdat_new[realdat_new$A==1,]

#Checking for extreme weights================================================================================

bb1 <- realdat_new %>% group_by(feedlot) %>% summarise(length(unique(A)))
dim(bb1[bb1[,2]!=1,]) # 64 cohorts contain more than 1 level of BVDV induction
temp <- as.matrix(bb1[bb1[,2]==1,][,1])
realdat_new_sub <- subset(realdat_new, feedlot %in% temp)
realdat_new_sub_0 <- realdat_new_sub[realdat_new_sub$A==0,]
realdat_new_sub_1 <- realdat_new_sub[realdat_new_sub$A==1,]

nrow(realdat_new_sub_0)
length(realdat_new_sub_0$icpw[realdat_new_sub_0$icpw==1])
length(realdat_new_sub_0$ps.fe[realdat_new_sub_0$ps.fe<0.01])
length(realdat_new_sub_0$ps.re[realdat_new_sub_0$ps.re<0.01])

nrow(realdat_new_sub_1)
length(realdat_new_sub_1$icpw[realdat_new_sub_1$icpw==1])
length(realdat_new_sub_1$ps.fe[realdat_new_sub_1$ps.fe>0.99])
length(realdat_new_sub_1$ps.re[realdat_new_sub_1$ps.re>0.99])

#================================================================================
# 
# dim(subdat_0)
# length(subdat_0$icpw[subdat_0$icpw==1])
# length(subdat_0$ps.fe[subdat_0$ps.fe<0.01])
# length(subdat_0$ps.re[subdat_0$ps.re<0.01])
# subdat_0[subdat_0$icpw>50,]
# 
# dim(subdat_1)
# length(subdat_1$icpw[subdat_1$icpw==1])
# length(subdat_1$ps.fe[subdat_1$ps.fe>0.99])
# length(subdat_1$ps.re[subdat_1$ps.re>0.99])
# subdat_1[subdat_1$icpw>20,]



