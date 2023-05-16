install.packages("Hmisc")
install.packages("epitools")
library(Hmisc)
library(epitools)
# use this code to bring in the dataset
vaccine<-read.table("PINKEYE.csv",sep=",",header=T,fill=T,  na.string=" ") 
# this code brings up the entire dataset- not very helpful
vaccine
# this code gives you descriptive information about the dataset
str(vaccine)
# this code lets you see that SEX is indeed a factor
is.factor(vaccine$SEX)
# this code lets you see that SWAB_NUMBER is not a factor
is.factor(vaccine$SWAB_NUMBER)
# You can calculate the mean of SWAB_NUMBER because its not a factor
mean(vaccine$SWAB_NUMBER)

# this code will only work if you have loaded Hmisc. This changes the variable names to lowercase
# this code also makes a new variable called swab where that is a factor
vaccine <- upData(vaccine, lowernames = TRUE, swab = factor(vaccine$SWAB_NUMBER))

str(vaccine)

mean(vaccine$swab)
# this code shows that can can't use swab as a number as its is a factor

ls()
# this code gives you a list of the dataset you have in R
# if you have more than one dataset you need to tell iR which to use
# one way is to always include the datsetname$variable name - this is slow
# note the head command returns the 1st n observations
head(vaccine$SEX, n = 15)

help(head)
# another way is to "attach" the dataset you want to use- this is not recommended 
search()
attach(vaccine)
search()
# notice that vaccine is now 1st the the search path 

detach(vaccine)
# this removes the attachment 
# another way to use specify the dataset to use is to use "with" 


with (vaccine ,(head(sex, n = 15)))

# lets start with descriptive continuous and binary data 
# no hypothesis testing yet

with (vaccine, (table(sex))) 
# create a table for the variable sex

help(table)

with (vaccine,(table(sex, producer)))  
# create a cross sectional table


mean(vaccine$weight) 
# should return a problem

mean(vaccine$weight, na.rm = TRUE)  
# will fix the problem

range(vaccine$weight, na.rm=TRUE)

with(vaccine, quantile(weight, probs = seq(0, 1, by = 0.1), na.rm=TRUE))

summary(vaccine$weight)
# this code will get you summary data for individual variables

names(vaccine)
# this is handy if you forget the names of your variables 


# its also possible to get summary tables if you have Hmisc loaded

summary(weight~vac_code, data = vaccine, methods= "responce")
summary(weight~vac_code, data = vaccine, methods= "responce", fun = function(x) {c(Median = median(x), Min = min(x), Max = max(x))})

summary(weight~vac_code, data = vaccine, method = "cross",fun = function(x) {
  c(Mean = mean(x), SD = sd(x), Median = median(x), Min = min(x),Max = max(x))})

summary(vac_code ~weight +sex+ producer, data = vaccine, method = "reverse", overall = TRUE)
# this code will provide you with a table that describe the basic descriptive information that you would need for most papers
print(summary(vac_code ~weight + sex + producer, overall = TRUE,
                      method = "reverse", data = vaccine), digits = 3, npct = "both",
      pctdig = 2, exclude1 = FALSE, long = TRUE, prmsd = TRUE) 

print(summary(pinkeye ~weight + sex + producer, overall = TRUE,
                      method = "reverse", data = vaccine), digits = 3, npct = "both",
      pctdig = 2, exclude1 = FALSE, long = TRUE, prmsd = TRUE) 

print(summary(pinkeye ~vac_code +weight+ sex, overall = TRUE,
                      method = "reverse", data = vaccine), digits = 3, npct = "both",
      pctdig = 2, exclude1 = FALSE, long = TRUE, test=TRUE, prmsd = TRUE) 

# What if you want to name some variables so they are more useful

table(vaccine$pinkeye, vaccine$vac_code)

vaccine$pinkeye<-factor(vaccine$pinkeye)
levels(vaccine$pinkeye)<-c("non-clinical", "clinical")
vaccine$vac_code<-factor(vaccine$vac_code)
levels(vaccine$vac_code)<-c("both","implant_only", "injectable_only", "control")

table(vaccine$pinkeye, vaccine$vac_code)

# testing some exploratory hypothesis testing for continuous data first for univariate analysis-
# This code will provide you with a linear model of weight compared to vac_code
fit = (lm (weight ~ vac_code, data=vaccine)) 
# this code asks for the analysis of variance
anova(fit) 

# this code asks for the analysis of variance
fit = (lm (weight ~ sex, data=vaccine)) 
anova(fit)

# this code asks for the analysis of variance
fit = (lm (weight ~ pinkeye, data=vaccine)) 
anova(fit)

# this code asks for the analysis of variance

fit = (lm (weight ~ producer, data=vaccine)) 
anova(fit)

# this code asks for the analysis of variance
fit = (lm (weight ~ location, data=vaccine)) 
anova(fit)

# now lets do exploratory analysis with categorical data : univariate 
# this code will get you table 
with (vaccine, table(sex,vac_code)) 
# this code will get you table 
with (vaccine, chisq.test(table, sex,vac_code))
with (vaccine, table(vac_code, pinkeye))
with (vaccine, chisq.test(vac_code, pinkeye)) 

# This will give you an over chi-square for the pinkeye
across the vaccine groups


# but what if you wanted to look at the analysis using only one producer
PROD1<-subset(vaccine,producer==1)
PROD2<-subset(vaccine,producer==2)
PROD3<-subset(vaccine,producer==3)
ls()

with (PROD1, table(vac_code, pinkeye))
with (PROD1, chisq.test(vac_code, pinkeye))
with (PROD2, table(vac_code, pinkeye))
with (PROD2, chisq.test(vac_code, pinkeye))
with (PROD3, table(vac_code, pinkeye))
with (PROD3, chisq.test(vac_code, pinkeye))

# this code asks for the analysis of variance
fit<-with (vaccine, lm(weight~sex + vac_code))
anova(fit)
# adding producer as a continuous variable
fit<-(lm (weight ~ sex + pinkeye + producer, data=vaccine)) 
anova(fit)

summary(weight~sex + pinkeye + producer, data = vaccine, method = "cross",fun = function(x) {
  c(Mean = mean(x))})


summary(fit)
summary(vaccine$producer)
vaccine <- upData(vaccine, lowernames = TRUE, treated = factor(vaccine$pinkeye), owner=factor(vaccine$producer))

fit<-(lm (weight ~ sex + treated + owner, data=vaccine)) 
summary(vaccine$owner)
anova(fit)
summary(fit)
plot(fitted(fit),residuals(fit), xlab="Fitted", ylab="residuals")

heifers<-subset(vaccine,sex==h)
bulls<-subset(vaccine,producer==2)
steers<-subset(vaccine,producer==3)
fit<-(lm(weight~sex, data=vaccine))
summary(fit)
model.matrix(fit)
#vaccine$sex<-relevel(vaccine$sex, ref="h")

vaccine$pinkeye<-factor(vaccine$pinkeye)
levels(vaccine$pinkeye)<-c("non-clinical", "clinical")
fit<-(lm(weight~pinkeye, data=vaccine))
summary(fit)
table(vaccine$pinkeye, vaccine$vac_code)
vaccine$vac_code<-factor(vaccine$vac_code)
levels(vaccine$vac_code)<-c("both","implant_only", "injectable_only", "control")

help(package=epitools)
help(epitab)
RRtab<-epitable(vaccine$vac_code, vaccine$pinkeye)
epitab(RRtab, method="riskratio", verbose = TRUE)
epitab(RRtab, method="oddsratio",  conf.level = 0.95, verbose = F)



dat.logistic= glm(pinkeye~vac_code,data=vaccine,family=binomial(logit))
summary(dat.logistic)
exp(0.21822)


