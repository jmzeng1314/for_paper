rm(list = ls())
load("database/icgc_tumor_normal_rawcount.Rdata")
library(tinyarray)
library(tidyverse)
exp[1:4,1:4]
identical(names(Group),colnames(exp))
exp_n=exp[,Group=='normal'];dim(exp_n)
exp_t = exp[,Group == "tumor"];dim(exp_t)

patient = str_split(colnames(exp_n),'-',simplify = T)[,1] #有配对样本的病人id
table((str_split(colnames(exp_t),'-',simplify = T)[,1]) %in% patient) # 看看肿瘤样本有重复吗

exp_t = exp_t[,(str_split(colnames(exp_t),'-',simplify = T)[,1]) %in% patient]
exp_t = exp_t[,str_sort(colnames(exp_t))]
exp_t = exp_t[,!duplicated(str_split(colnames(exp_t),'-',simplify = T)[,1])] #排序去重，这样样本有01A和01B可选时，就会留下01A。
dim(exp_t)

#让正常样本顺序与肿瘤样本顺序一致,这样才方便写配对向量
tmp = match((str_split(colnames(exp_t),'-',simplify = T)[,1]),(str_split(colnames(exp_n),'-',simplify = T)[,1]));head(tmp)

exp_n = exp_n[,tmp]
identical((str_split(colnames(exp_t),'-',simplify = T)[,1]),(str_split(colnames(exp_n),'-',simplify = T)[,1]))
table(str_split(colnames(exp_t),'-',simplify = T)[,2])
table(str_split(colnames(exp_n),'-',simplify = T)[,2])
exp_small = cbind(exp_n,exp_t)
Group_small = rep(c("normal","tumor"),each = ncol(exp_n))
Group_small = factor(Group_small,levels = c("normal","tumor"))
pairinfo = rep(1:ncol(exp_n),times = 2)

table(rownames(cl_tumor)%in%patient)
clinical_small=cl_tumor[rownames(cl_tumor)%in%patient,];dim(clinical_small)
save(exp_small,Group_small,clinical_small,file = 'database/icgc_pair_exp_clinical.Rdata')
