rm(list = ls())
load('database/tcga_pair_exp_clinical.Rdata')
clinical=clinical_small
meta = clinical[,c(
  'bcr_patient_barcode',
  'vital_status',
  'days_to_death',
  'days_to_last_followup',
  'race_list',
  'days_to_birth',
  'gender' ,
  'stage_event',
  'relative_family_cancer_history',
  'history_hepato_carcinoma_risk_factors',
  'neoplasm_histologic_grade',
  'fetoprotein_outcome_value',
  'platelet_result_count',
  'prothrombin_time_result_value',
  'albumin_result_specified_value',
  'fibrosis_ishak_score',
  'adjacent_hepatic_tissue_inflammation_extent_type',
  'viral_hepatitis_serologies',
  'new_tumor_events'
)]
dim(meta)
rownames(meta) <- meta$bcr_patient_barcode
meta[1:4,1:4]
colnames(meta)=c('ID','event','death','last_followup','race','age','gender','stage',
                 'family_history','causes','inflammation_grade','AFP','PLT',
                 'PT','ALB','fibrosis_score','adjacent_inflammation',
                 'viral_hepatitis_serologies','new_tumor_events')
#空着的值改为NA
meta[meta==""]=NA
# 1.1由随访时间和死亡时间计算生存时间(月)
table(meta$event)
meta$time = ifelse(meta$event=="Alive",
                   meta$last_followup,
                   meta$death)
meta$time = as.numeric(meta$time)/30
#1.2 根据生死定义event，活着是0，死的是1
table(meta$event)
# 1.3 年龄和年龄分组
meta$age=ceiling(abs(as.numeric(meta$age))/365)

meta$age_group=ifelse(meta$age>65,'older','younger')
table(meta$age_group)
# 1.4 stage 
library(stringr) 
head(meta$stage)
a=str_split(meta$stage,'T',simplify=T)
a=str_split(a[,2],'',simplify=T)
a=a[,1]
meta$T_stage=a
a=str_split(meta$stage,'N',simplify=T)
a=str_split(a[,2],'',simplify=T)
a=a[,1]
meta$N_stage=a
a=str_split(meta$stage,'M',simplify=T)
a=str_split(a[,2],'',simplify=T)
a=a[,1]
meta$M_stage=a
table(meta$M_stage)
meta$T_stage=ifelse(meta$T_stage=='1',1,
                    ifelse(meta$T_stage=='2',2,
                           ifelse(meta$T_stage=='3',3,
                                  ifelse(meta$T_stage=='4',4,
                                         NA))))
meta$N_stage=ifelse(meta$N_stage=='0',0,
                    ifelse(meta$N_stage=='1',1,NA
                    ))

meta$M_stage=ifelse(meta$M_stage=='0',0,
                    ifelse(meta$M_stage=='1',1,NA
                    ))

a = str_extract_all(meta$stage,"I|V");head(a)
b = sapply(a,paste,collapse = "");head(b)
meta$stage = b

str(meta)
colnames(meta)
meta=meta[,c(1,2,6,7,8,9,11,20,21)]
str(meta)
meta$stage=factor(meta$stage,levels = c('I','II','III','IV'))
meta$event=factor(meta$event,levels = c('Alive','Dead'))
meta$age_group=factor(meta$age_group,levels = c('younger','older'))
meta$family_history=factor(meta$family_history,levels = c('NO','YES'))
meta$gender=ifelse(meta$gender=='MALE','Male','Female')
meta$gender=factor(meta$gender,levels = c('Male','Female'))
table(meta$inflammation_grade)
save(exp_small,meta,file = 'Rdata/tcga/tcga_pair_exp_meta_zhenglihou.Rdata')



load('database/icgc_tumor_normal_rawcount.Rdata')
colnames(exp)[1:10]
exp_normal=exp[,str_detect(colnames(exp),'normal')]
rownames(cl_tumor)
a=str_split(colnames(exp_normal),'-',simplify = T)[,1]
cl=cl_tumor[rownames(cl_tumor)%in%a,]
str(cl)
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
colnames(meta)
cl$family_history=rep(NA,nrow(cl))
cl$inflammation_grade=rep(NA,nrow(cl))
colnames(cl)
cl$group=rep('ICGC',nrow(cl))
meta$group=rep('TCGA',nrow(meta))
a=c('age','age_group','gender','family_history','inflammation_grade','stage','time','event','group')
cl=cl[,match(a,colnames(cl))]
meta=meta[,match(a,colnames(meta))]
identical(colnames(meta),colnames(cl))
dat=rbind(meta,cl)
dat$group=factor(dat$group,levels = c('TCGA','ICGC'))
save(dat,file = 'database/pair_icgc_tcga_clincal_zhengli.Rdata')
table=table1(~ gender + age+age_group +  family_history+
               inflammation_grade+stage +time+event|group, data=dat, overall="Total")
