rm(list = ls())
load('Rdata/tcga/pre/clinicaldata_juzhenshuju.Rdata')
colnames(meta)
meta=meta[,c(2,6,7,8,9,11,20,21,26,27)]
str(meta)
meta$stage=factor(meta$stage,levels = c('I','II','III','IV'))
meta$event=ifelse(meta$event==0,'Alive','Dead')
meta$event=factor(meta$event,levels = c('Alive','Dead'))
meta$age_group=factor(meta$age_group,levels = c('younger','older'))
meta$family_history=factor(meta$family_history,levels = c('NO','YES'))
meta$gender=ifelse(meta$gender=='MALE','Male','Female')
meta$gender=factor(meta$gender,levels = c('Male','Female'))
table(meta$inflammation_grade)
meta_tcga=meta
load('Rdata/icgc/pre/linchuangzuiquan_icgc.Rdata')
colnames(meta)
cl=meta
cl$age_group=ifelse(cl$age<65,'younger','older')
cl$age_group=factor(cl$age_group,levels=c('younger','older'))
cl$stage=ifelse(cl$stage==1,'I',
                ifelse(cl$stage==2,'II',
                       ifelse(cl$stage==3,'III','IV')))
cl$stage=factor(cl$stage,levels = c('I','II','III','IV'))
cl$event=ifelse(cl$event==0,'Alive','Dead')
cl$event=factor(cl$event,levels = c('Alive','Dead'))
cl$gender=ifelse(cl$gender=='male','Male','Female')
class(cl$gender)
cl$gender=factor(cl$gender,levels = c('Male','Female'))
colnames(cl)
colnames(meta_tcga)
cl$family_history=rep(NA,nrow(cl))
cl$inflammation_grade=rep(NA,nrow(cl))
colnames(cl)
cl$group=rep('ICGC',nrow(cl))
meta_tcga$group=rep('TCGA',nrow(meta_tcga))
colnames(cl)
cl=cl[,-6]
colnames(meta_tcga)
a=c('age','age_group','gender','family_history','inflammation_grade','stage','time','event','riskscore','ri','group')
cl=cl[,match(a,colnames(cl))]
meta_tcga=meta_tcga[,match(a,colnames(meta_tcga))]
identical(colnames(meta_tcga),colnames(cl))
dat=rbind(meta_tcga,cl)
dat$group=factor(dat$group,levels = c('TCGA','ICGC'))
save(dat,file = 'Rdata/tumor_icgc_tcga_clincal_zhengli.Rdata')
table=table1(~ gender + age+age_group +  family_history+
               inflammation_grade+stage +time+riskscore+ri|group*event, data=dat, overall="Total")
