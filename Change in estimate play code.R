#	Example on page 146-152 of Logistic Regression by Kleinbaum

#Load data===============
evans<-read.table("evans.dat",header=F,col.names=c("ID","CHD","CAT","AGE",
                                                   "CHL","SMK","ECG","DBP","SBP","HPT","CH","CC"))
#	Fit a logistic regression model without interaction========
evans.fit1<-glm(CHD~CAT+AGE+CHL+ECG+SMK+HPT,family=binomial,data=evans)

#	Print summary information============
summary(evans.fit1)		#printout in page 146 and 148
est.b2<-round(cbind(evans.fit1$coef,confint.default(evans.fit1)),2)[2,]
est.b2  	# beta estimate in page 147

#	Odds ratio estimates   
or=exp(est.b2)
or 
round(or,2)	# OR estimate in page 147

#residual deviance model without interaction
#P values ========
# #wald test page 148
# Z<-evans.fit1$coef[[2]]/sqrt(diag(vcov(evans.fit1)))[[2]]
# Z
# Z2<-Z^2
# Z2
# 
# #p-value (assume two tailed test)
# 1-pchisq(Z2,1)
# 
# #p-value one tailed
# (1-pchisq(Z2,1))/2                                                    

#Estimates : logOR and OR with 95% CI interval=======
# beta estimates and 95% confidence limits page 149 
est.b<-round(cbind(evans.fit1$coef,confint.default(evans.fit1)),2)
est.b		# beta estimate in page 147

#OR estimates and 95% confidence limits
round(exp(est.b),2)     