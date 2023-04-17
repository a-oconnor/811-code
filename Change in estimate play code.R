#	Example on page 146-152

#	Load data
evans<-read.table("evans.dat",header=F,col.names=c("ID","CHD","CAT","AGE",
                                                   "CHL","SMK","ECG","DBP","SBP","HPT","CH","CC"))

#	Fit a logistic regression model without interaction
evans.fit1<-glm(CHD~CAT+AGE+CHL+ECG+SMK+HPT,family=binomial,data=evans)
#	Print summary information
summary(evans.fit1)		#printout in page 146 and 148


est.b2<-round(cbind(evans.fit1$coef,confint.default(evans.fit1)),2)[2,]
est.b2  	# beta estimate in page 147

#	Odds ratio estimates   
or=exp(est.b2)
or 
round(or,2)	# OR estimate in page 147

#residual deviance model without interaction

#wald test page 148
Z<-evans.fit1$coef[[2]]/sqrt(diag(vcov(evans.fit1)))[[2]]
Z
Z2<-Z^2
Z2

#p-value (assume two tailed test)
1-pchisq(Z2,1)

#p-value one tailed
(1-pchisq(Z2,1))/2                                                    


# beta estimates and 95% confidence limits page 149 
est.b<-round(cbind(evans.fit1$coef,confint.default(evans.fit1)),2)
est.b		# beta estimate in page 147

#OR estimates and 95% confidence limits
round(exp(est.b),2)                                                     

#	Fit a logistic regression model with interaction pag 147
evans.fit2=glm(CHD~CAT+AGE+CHL+ECG+SMK+HPT+CH+CC,family=binomial,data=evans)
#	Print summary information
summary(evans.fit2)		#printout in page 147

#residual deviance
evans.fit2$deviance

#	LRT
lr<-evans.fit1$deviance-evans.fit2$deviance
lr
1-pchisq(lr,2)		#LRT in page 150

#	ROR for CAT=1 vs. CAT=0 for HPT=0,1 and CHL=200,220,240 page 151

cc<-matrix(c(0,1,0,0,0,0,0,0,200,	#HPT=0, CHL=200
             0,1,0,0,0,0,0,0,220,	#HPT=0, CHL=220
             0,1,0,0,0,0,0,0,240,	#HPT=0, CHL=240
             0,1,0,0,0,0,0,1,200,	#HPT=1, CHL=200
             0,1,0,0,0,0,0,1,220,	#HPT=1, CHL=220
             0,1,0,0,0,0,0,1,240),6,9,byrow=T)	#HPT=1, CHL=240

co2<-round(evans.fit2$coef,4)
est<-cc%*%co2
round(exp(est),2)


round(est[5],4)						#l.hat
vcov(evans.fit2)				#variance covariance matrix
v.est<-cc%*%vcov(evans.fit2)%*%t(cc)	
round(v.est[5,5],4)					#v.hat(l.hat)

se<-sqrt(diag(cc%*%vcov(evans.fit2)%*%t(cc)))
LL<-est-qnorm(.975)*se
UL<-est+qnorm(.975)*se
ROR<-data.frame(Estimate=exp(est),LL=exp(LL),UL=exp(UL))
rownames(ROR)<-c("HPT0CHL200","HPT0CHL220","HPT0CHL240","HPT1CHL200","HPT1CHL220","HPT1CHL240")
round(ROR,2)		#OR and 95% CI in page 151 and 152

