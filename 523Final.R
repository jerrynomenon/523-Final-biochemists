library(tidyverse) 
library(stargazer)
library(xtable)

file<-"biochemists.csv"
data<-read.csv(file)

art <- data$art
fem <- data$fem
mar <- data$mar
kid5 <- data$kid5
phd <- data$phd
ment <- data$ment
anypub <- data$anypub

# 1
# all three are factors with 2 lvls - use 2x2 table

#poisson with table
table(fem,mar,anypub)
as.data.frame(table(fem,mar,anypub))
fem.count<-as.data.frame(table(fem,mar,anypub))$fem
mar.count<-as.data.frame(table(fem,mar,anypub))$mar
anypub.count<-as.data.frame(table(fem,mar,anypub))$anypub
count<-as.data.frame(table(fem,mar,anypub))$Freq

fit.fma<-glm(count~fem.count*mar.count*anypub.count,data=data,family="poisson")
summary(fit.fma)

# high p value for f:m:a
fit.fm.ma.fa<-update(fit.fma, ~.-fem.count:mar.count:anypub.count)
summary(fit.fm.ma.fa)

fit.fm.ma<-update(fit.fm.ma.fa, ~.- fem.count:anypub.count)
summary(fit.fm.ma)

fit.fm.fa<-update(fit.fm.ma.fa, ~.- mar.count:anypub.count)
summary(fit.fm.ma)

# ma and fa both have high p val

fit.fm.a<-update(fit.fm.fa, ~. - fem.count:anypub.count)
summary(fit.fm.a)

fit.f.m.a<-update(fit.fm.a, ~. - fem.count:mar.count)
summary(fit.f.m.a)

data.frame(as.data.frame(table(fem,mar,anypub)), fma = fitted(fit.fma),fm.ma.fa = fitted(fit.fm.ma.fa),fm.ma = fitted(fit.fm.ma),fm.fa = fitted(fit.fm.fa), fm.a = fitted(fit.fm.a), f.m.a = fitted(fit.f.m.a))

# compare models against saturated model using goodness of fit
modellist = list(fit.f.m.a, fit.fm.a, fit.fm.fa, fit.fm.ma, fit.fm.ma.fa, fit.fma)

data.frame(Model=c("f.m.a", "fm.a", "fm.fa", "fm.ma", "fm.ma.fa", "fma"), Dev=round(unlist(lapply(modellist, deviance)),4),
           X2=round(unlist(lapply(modellist, function(x){ sum( residuals(x,type="pearson")^2)})),4),
           AIC=round(unlist(lapply(modellist, AIC)),4),
           BIC=round(unlist(lapply(modellist, AIC,k=log(915))),4),
           Df = unlist(lapply(modellist,function(x){ x$df.residual})))
# evidence against f.m.a and fm.ma
#X2 test
anova(fit.fm.a, fit.fm.fa, fit.fm.ma.fa, fit.fma,test="Chisq")
anova(fit.fm.a, fit.fm.fa, fit.fm.ma.fa, fit.fma,test="LRT")
# fit.fm.a is most suitable model 
# conf int
confint(fit.fm.a)
exp(confint(fit.fm.a))
plot(fit.fm.a)

#2
pairs(cbind(art, phd, ment))

# art and phd seems to have no significant correlation, though phd and ment seem to do
par(mfrow=c(1,1))
boxplot(art~fem)
# men tends to publish more articles.
boxplot(art~mar)
# Married people to publish more articles, makes sense as married person tends to be older, thus more experienced and renown.
boxplot(art~kid5)
# People with less kids tend to publish more articles.

# 

fit.f.m.k.p.mt<-glm(art ~ fem + mar + kid5 + phd + ment, family = poisson, data = data)
summary(fit.f.m.k.p.mt)

# add interaction between phd and mentor
fit.f.m.k.pmt<-glm(art ~ fem + mar + kid5 + phd * ment, family = poisson, data = data)
summary(fit.f.m.k.pmt)

# or take out phd
fit.f.m.k.mt<-update(fit.f.m.k.p.mt, ~. - phd)
summary(fit.f.m.k.mt)

# take out mar
fit.f.k.mt<-update(fit.f.m.k.mt, ~. - mar)
summary(fit.f.k.mt)

pairs(cbind(kid5, ment, fem))
# kid5 and mentor have some correlation
fit.f.m.kmt.pmt<-update(fit.f.m.k.pmt, ~. + kid5:ment)
summary(fit.f.m.kmt.pmt)

modellist2 = list(fit.f.k.mt, fit.f.m.k.mt, fit.f.m.k.p.mt, fit.f.m.k.pmt, fit.f.m.kmt.pmt)
data.frame(Model=c("f.k.mt", "f.m.k.mt", "f.m.k.p.mt", "f.m.k.pmt", "fit.f.m.kmt.pmt"), Dev=round(unlist(lapply(modellist, deviance)),4),
      X2=round(unlist(lapply(modellist2, function(x){ sum( residuals(x,type="pearson")^2)})),4),
      AIC=round(unlist(lapply(modellist2, AIC)),4),
      BIC=round(unlist(lapply(modellist2, AIC,k=log(915))),4),
      Df = unlist(lapply(modellist2,function(x){ x$df.residual})))

anova(fit.f.k.mt, fit.f.m.k.mt, fit.f.m.k.p.mt, fit.f.m.k.pmt, fit.f.m.kmt.pmt,test="LRT")
# significant enough correlation between phd and ment, keep


X2 =sum(residuals(fit.f.m.kmt.pmt,type="pearson")^2)
cat(c("X2 = ",round(X2,3)))
cat(c("Rejection value = ",round(qchisq(0.95,914),3)))
phi = X2/(907)
cat(c("Phi=",round(phi,3)))
# phi is 1.8, let's see how nb performs
fit.f.m.kmt.pmt.nb<-glm.nb(art ~ fem + mar + kid5 + phd + ment + phd:ment + kid5:ment, data = data)
summary(fit.f.m.kmt.pmt.nb)
fit.f.m.kmt.pmt.quasi<-glm(art ~ fem + mar + kid5 + phd + ment + phd:ment + kid5:ment,family=quasi(link="log",variance="mu"), data = data)
summary(fit.f.m.kmt.pmt.quasi)
fit.f.m.kmt.pmt.quasi2<-glm(art ~ fem + mar + kid5 + phd + ment + phd:ment + kid5:ment,family=quasi(link="log",variance="mu^2"), data = data)
summary(fit.f.m.kmt.pmt.quasi2)
# neither quasi nor nb performs well


