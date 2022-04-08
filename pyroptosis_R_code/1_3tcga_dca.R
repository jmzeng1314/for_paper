rm(list=ls())
load('Rdata/tcga/pre/clinicaldata_juzhenshuju.Rdata')

library(regplot)
library(survival)
identical(colnames(exp_dat),rownames(meta))
meta=meta[match(colnames(exp_dat),rownames(meta)),]
identical(colnames(exp_dat),rownames(meta))
dat=cbind(meta,t(exp_dat))
colnames(dat)

dat$gender=ifelse(dat$gender=='MALE',0,1)
dat$stage=ifelse(dat$stage=='I',1,
                 ifelse(dat$stage=='II',2,
                        ifelse(dat$stage=='III',3,4)))
dat$family_history=ifelse(dat$family_history=='NO',0,1)                 
dat$inflammation_grade=ifelse(dat$inflammation_grade=='G1',1,
                              ifelse(dat$inflammation_grade=='G2',2,
                                     ifelse(dat$inflammation_grade=='G3',3,4)))
dat$adjacent_inflammation=ifelse(dat$adjacent_inflammation=='None',0,
                                 ifelse(dat$adjacent_inflammation=='Mild',1,2)
)
library(rms)
dd<-datadist(dat)
options(datadist="dd")
library(ggDCA)
library(rms)
colnames(dat)
riskscore <- cph(formula = Surv(time, event) ~ riskscore,
                 data=dat,x=T,y=T,surv = T)
stage <- cph(formula = Surv(time, event) ~ stage,
             data=dat,x=T,y=T,surv = T)
inflammation_grade <- cph(formula = Surv(time, event) ~ inflammation_grade,
                          data=dat,x=T,y=T,surv = T)
#adjacent_inflammation <- cph(formula = Surv(time, event) ~ adjacent_inflammation,
#data=dat,x=T,y=T,surv = T)
#T_stage <- cph(formula = Surv(time, event) ~ T_stage,
               #data=dat,x=T,y=T,surv = T)
age <- cph(formula = Surv(time, event) ~ age,
           data=dat,x=T,y=T,surv = T)
gender <- cph(formula = Surv(time, event) ~ gender,
              data=dat,x=T,y=T,surv = T)

#family_history <- cph(formula = Surv(time, event) ~ family_history,
                      #data=dat,x=T,y=T,surv = T)

data  <- dca(riskscore,
             stage,inflammation_grade,age,gender)
ggplot(data, 
       linetype =T)   #线粗
ggsave(file='plot/tcga/pre/dea_tcga.png',height = 5,width = 5)
