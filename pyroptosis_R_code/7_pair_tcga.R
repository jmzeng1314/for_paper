rm(list = ls())
load("database/TCGA-LIHC_gdc.Rdata")
library(tinyarray)
library(tidyverse)
exp[1:4,1:4]
exp_n=exp[,Group=='normal'];dim(exp_n)
exp_t = exp[,Group == "tumor"];dim(exp_t)

patient = str_sub(colnames(exp_n),1,12) #有配对样本的病人id
table(str_sub(colnames(exp_t),1,12) %in% patient) # 看看肿瘤样本有重复吗

exp_t = exp_t[,str_sub(colnames(exp_t),1,12) %in% patient]
exp_t = exp_t[,str_sort(colnames(exp_t))]
exp_t = exp_t[,!duplicated(str_sub(colnames(exp_t),1,12) )] #排序去重，这样样本有01A和01B可选时，就会留下01A。
dim(exp_t)

#让正常样本顺序与肿瘤样本顺序一致,这样才方便写配对向量
tmp = match(str_sub(colnames(exp_t),1,12),str_sub(colnames(exp_n),1,12));head(tmp)

exp_n = exp_n[,tmp]
identical(str_sub(colnames(exp_t),1,12),str_sub(colnames(exp_n),1,12))
table(str_sub(colnames(exp_t),14,15))
table(str_sub(colnames(exp_n),14,15))
exp_small = cbind(exp_n,exp_t)
Group_small = rep(c("normal","tumor"),each = ncol(exp_n))
Group_small = factor(Group_small,levels = c("normal","tumor"))
pairinfo = rep(1:ncol(exp_n),times = 2)

table(clinical$bcr_patient_barcode%in%patient)
clinical_small=clinical[clinical$bcr_patient_barcode%in%patient,];dim(clinical_small)
save(exp_small,Group_small,clinical_small,file = 'database/tcga_pair_exp_clinical.Rdata')
