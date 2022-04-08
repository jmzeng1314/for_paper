rm(list = ls())
load("database/TCGA-LIHC_gdc.Rdata")
library(tinyarray)
library(tidyverse)
exp[1:4,1:4]
exp = trans_exp(exp,mrna_only = T)
exp[1:4,1:4]
Group=ifelse(str_sub(colnames(exp),14,15)=='01','tumor','normal')
Group=factor(Group,levels = c('normal','tumor'))
exp=exp[,Group=='tumor']
exp=log2(cpm(exp)+1)
dim(exp)
#与细胞焦亡相关基因取交集
load('database/pyroptosis.Rdata')
table(rownames(exp)%in%pyroptosis)
pyroptosis[!(pyroptosis%in%rownames(exp))]
exp=exp[rownames(exp)%in%pyroptosis,]
dim(exp)
save(exp,clinical,Group,file='Rdata/tcga/pre/pyroptosis_tumor.Rdata')
dim(exp)
exprSet = exp
dim(exprSet)
tmp = data.frame(colnames(clinical))
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
# 以病人为中心，表达矩阵按病人ID去重复
k = !duplicated(str_sub(colnames(exprSet),1,12));table(k)

exprSet = exprSet[,k]
#调整meta的ID顺序与exprSet列名一致
meta=meta[match(str_sub(colnames(exprSet),1,12),meta$ID),]
identical(meta$ID,str_sub(colnames(exprSet),1,12))
# 1.1由随访时间和死亡时间计算生存时间(月)
table(meta$event)
meta$time = ifelse(meta$event=="Alive",
                   meta$last_followup,
                   meta$death)
meta$time = as.numeric(meta$time)/30
#1.2 根据生死定义event，活着是0，死的是1
meta$event=ifelse(meta$event=='Alive',
                  0,
                  1)
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

# 去掉生存信息不全或者生存时间小于0.1(月)的病人，样本纳排标准不唯一，且差别很大
k1 = meta$time>=0.1;table(k1)
k2 = !(is.na(meta$time)|is.na(meta$event));table(k2)

exprSet = exprSet[,k1&k2]
meta = meta[k1&k2,]
table(meta$stage)
meta=meta[meta$stage!='',]
summary(meta$time)
exprSet=exprSet[,str_sub(colnames(exprSet),1,12)%in%rownames(meta)]
identical(str_sub(colnames(exprSet),1,12),rownames(meta))

save(meta,exprSet,file = "Rdata/tcga/pre/TCGA-LIHC_sur_model.Rdata")

rm(list = ls())
load("Rdata/tcga/pre/TCGA-LIHC_sur_model.Rdata")
library(survival)
library(survminer)
logrankfile = "Rdata/tcga/pre/log_rank_p.Rdata"


if(!file.exists(logrankfile)){
  log_rank_p <- apply(exprSet , 1 , function(gene){
    # gene=exprSet[1,]
    meta$group=ifelse(gene>median(gene),'high','low')  
    data.survdiff=survdiff(Surv(time, event)~group,data=meta)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    return(p.val)
  })
  log_rank_p=sort(log_rank_p)
  save(log_rank_p,file = logrankfile)
}
load(logrankfile)
table(log_rank_p<0.01) 
table(log_rank_p<0.05) 

coxfile ="Rdata/tcga/pre/unti_cox.Rdata"
if(!file.exists(coxfile)){
  cox_results <-apply(exprSet , 1 , function(gene){
    #gene= exprSet[1,]
    meta$gene = gene
    #可直接使用连续型变量
    m = coxph(Surv(time, event) ~ gene, data =  meta)
    #也可使用二分类变量
    #meta$group=ifelse(gene>median(gene),'high','low') 
    #m=coxph(Surv(time, event) ~ group, data =  meta)
    beta <- coef(m)
    se <- sqrt(diag(vcov(m)))
    HR <- exp(beta)
    HRse <- HR * se
    
    #summary(m)
    tmp <- round(cbind(coef = beta, 
                       se = se, z = beta/se, 
                       p = 1 - pchisq((beta/se)^2, 1),
                       HR = HR, HRse = HRse,
                       HRz = (HR - 1) / HRse, 
                       HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                       HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                       HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
    
    return(tmp['gene',]) 
    #return(tmp['grouplow',])#二分类变量
  })
  cox_results=as.data.frame(t(cox_results))
  save(cox_results,file = coxfile)
}
load(coxfile)
table(cox_results$p<0.01)
table(cox_results$p<0.2)
lr = names(log_rank_p)[log_rank_p<0.05];length(lr)
cox = rownames(cox_results)[cox_results$p<0.2];length(cox)
length(intersect(lr,cox))


exp_cox_gene=exprSet[rownames(exprSet)%in%cox,]
dim(exp_cox_gene)
save(cox,exp_cox_gene,meta,file = 'Rdata/tcga/pre/cox_gene_data.Rdata')
lr_cox_gene=union(lr,cox);length(lr_cox_gene)
exp_lr_cox_gene=exprSet[rownames(exprSet)%in%lr_cox_gene,]
table(rownames(exprSet)%in%lr_cox_gene)
dim(exp_lr_cox_gene)
save(lr_cox_gene,exp_lr_cox_gene,meta,file = 'Rdata/tcga/pre/lr_cox_gene_data.Rdata')
a=cox_results
a=a[rownames(a)%in%cox,]
write.csv(a,file = 'export/tcga/pre/unicox_tcga.csv')
