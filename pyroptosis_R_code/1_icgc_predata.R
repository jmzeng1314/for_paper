rm(list = ls())
# pd=read.csv('import/donor.LIRI-JP.tsv/donor.LIRI-JP.tsv',sep = '')
# exp=read.csv('import/exp_seq.LIRI-JP.tsv/exp_seq.LIRI-JP.tsv',sep='\t')
# sp=read.csv('import/specimen.LIRI-JP.tsv/specimen.LIRI-JP.tsv',sep='\t')
# save(pd,exp,sp,file = 'export/icgc_raw.data')
load('export/icgc_raw.data')
#安装R包
library(readr)
library(dplyr)
library(tidyverse)
library(tibble)

#（二）提取癌症的specimen_ID & donor_ID
#https://www.yuque.com/annananaa/geifhs/wkf6kz(注意参考资料中是用read_tsv进行读取文件的，这里读取的文件可以是未解压的文件)
disease='LIHC'
specimen=sp
specimen_ID = specimen[specimen$specimen_type == "Primary tumour - solid tissue",]$icgc_specimen_id
exp = exp[exp$icgc_specimen_id %in% specimen_ID,]
colnames(exp)
head(exp)
donor_ID = unique(exp$icgc_donor_id)
donor=pd
colnames(donor)
donor = donor[donor$icgc_donor_id %in% donor_ID,]

donor1 = donor[,c("icgc_donor_id","donor_sex","disease_status_last_followup",
                 "donor_vital_status","donor_diagnosis_icd10" ,
                 "donor_tumour_staging_system_at_diagnosis")]


colnames(donor1)=c("icgc_donor_id","gender","age",
                 "event","stage",
                 "time")

length(unique(exp$icgc_donor_id)) 
as.numeric(table(exp$icgc_donor_id)) 

# (三）表达矩阵整理
# icgc是长矩阵的形式储存信息，所以基因都保存到了一列，要把它们整理为常见的形式，无脑运行下面的
dim(donor1)
dim(exp)

colnames(exp)

exp = exp[,c("icgc_donor_id", "gene_id", "raw_read_count" )]

#以上三句代码可以用管道符号链接起来
if(T){
  exp = exp %>% group_by(icgc_donor_id, gene_id) %>% 
    summarise_all(max) %>% 
    pivot_wider(names_from = "icgc_donor_id", values_from = "raw_read_count") %>%
    summarise_all(function(x){ifelse(is.na(x),0,x)}) 
}
dim(exp)
dim(donor1)

exp = exp[!duplicated(exp$gene_id),]
exp = tibble::column_to_rownames(exp, "gene_id")

donor1=donor1[donor1$time!='yes',]
donor1=donor1[donor1$time!='no',]
dim(donor)
dim(donor1)
a=donor$icgc_donor_id[!donor$icgc_donor_id%in%donor1$icgc_donor_id]
a=donor[!donor$icgc_donor_id%in%donor1$icgc_donor_id,]
colnames(a)
a=a[,c(1,4,5,6,12,13)]
colnames(a)=c("icgc_donor_id","gender",
              "event","age","stage",
              "time")
a=a[,match(colnames(donor1),colnames(a))]
donor=rbind(donor1,a)
dim(donor)
exp=exp[,colnames(exp)%in%donor$icgc_donor_id]
dim(exp)
cl = donor[match(colnames(exp),donor$icgc_donor_id),]
identical(cl$icgc_donor_id,colnames(exp))
b=cl$icgc_donor_id
rownames(cl)=cl$icgc_donor_id
cl=cl[,-1]

cl$event = ifelse(cl$event == "alive", 0, 1)
str(cl)
cl$gender=factor(cl$gender,levels = c('male','female'))
cl$age=as.numeric(cl$age)
cl$stage=as.numeric(cl$stage)
cl$time=as.numeric(cl$time)
cl$time=cl$time/30
save(exp, cl, file = paste0("export/icgc_", disease, "_exp_cl.Rdata"))
