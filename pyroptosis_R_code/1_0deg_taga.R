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
Group_small = rep(c("non-tumor","tumor"),each = ncol(exp_n))
Group_small = factor(Group_small,levels = c("non-tumor","tumor"))
pairinfo = rep(1:ncol(exp_n),times = 2)

table(clinical$bcr_patient_barcode%in%patient)
clinical_small=clinical[clinical$bcr_patient_barcode%in%patient,];dim(clinical_small)
save(exp_small,Group_small,clinical_small,file = 'Rdata/tcga/pair_exp_clinical.Rdata')


load('Rdata/tcga/pair_exp_clinical.Rdata')
exp=exp_small
Group=Group_small
exp = trans_exp(exp,mrna_only = T)
exp[1:4,1:4]
a=exp
#与细胞焦亡相关基因取交集
load('database/pyroptosis.Rdata')
table(rownames(exp)%in%pyroptosis)
pyroptosis[!(pyroptosis%in%rownames(exp))]
exp=exp[rownames(exp)%in%pyroptosis,]
dim(exp)
save(exp,clinical,Group,file='Rdata/pyroptosis_tumor_normal.Rdata')
table(Group)
write.table(pyroptosis,file = 'export/51_genes.txt',quote =F,row.names = F)

library(DESeq2)
colData <- data.frame(row.names =colnames(exp), 
                      condition=Group)
if(!file.exists("Rdata/tcga/pyroptosis_DESeq.Rdata")){
  dds <- DESeqDataSetFromMatrix(
    countData = exp,
    colData = colData,
    design = ~ condition)
  dds <- DESeq(dds)
  save(dds,file = "Rdata/tcga/pyroptosis_DESeq.Rdata")
}
load(file ="Rdata/tcga/pyroptosis_DESeq.Rdata")
class(dds)
res <- results(dds, contrast = c("condition",rev(levels(Group))))
c("condition",rev(levels(Group)))

class(res)

DEG1 <- as.data.frame(res)
DEG1 <- DEG1[order(DEG1$pvalue),] 
DEG1 = na.omit(DEG1)
head(DEG1)


#添加change列标记基因上调下调
logFC_t = 0
pvalue_t = 0.05

k1 = (DEG1$padj< pvalue_t)&(DEG1$log2FoldChange < -logFC_t);table(k1)
k2 = (DEG1$padj < pvalue_t)&(DEG1$log2FoldChange > logFC_t);table(k2)
DEG1$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG1$change)
head(DEG1)

save(DEG1,Group_small,file = "Rdata/tcga/DEG1.Rdata")
write.csv(DEG1,file = 'export/tcga/deseq2_deg.csv')

#可视化
library(ggplot2)
library(tinyarray)
exp=a
dim(exp)
exp[1:4,1:4]
# cpm 去除文库大小的影响
dat = log2(cpm(exp)+1)
load('database/pyroptosis.Rdata')
table(rownames(dat)%in%pyroptosis)
pyroptosis[!(pyroptosis%in%rownames(dat))]
dat=dat[rownames(dat)%in%pyroptosis,]
dim(dat)
b=rownames(DEG1)[DEG1$change!='NOT']
table(Group)
Group=ifelse(str_sub(colnames(dat),14,15)=='01','tumor','non-tumor')
Group=factor(Group,levels = c('non-tumor','tumor'))
table(Group)
pca.plot = draw_pca(dat[rownames(dat)%in%b,],Group,addEllipses = F);pca.plot

save(pca.plot,file = "Rdata/tcga/tumor_normal_pcaplot.Rdata")
ggsave('plot/tcga_pair_pca.png',height = 5,width = 5)

cg1 = rownames(DEG1)[DEG1$change !="NOT"]

h1 = draw_heatmap(dat[cg1,],Group,n_cutoff = 2,cluster_cols = F,show_rownames = T)
h2 = draw_heatmap(dat[cg1,],Group,n_cutoff = 2,cluster_cols = F,show_rownames = T,annotation_legend = T,legend = T)
dev.off()
h2
ggsave('plot/tcga_pair_heatmap.png',height = 5,width = 5)
table(DEG1$change)
k=DEG1
v1 = draw_volcano(k,pkg = 1,logFC_cutoff = logFC_t)
table(DEG1$change)
#维恩图拼图，终于搞定
library(patchwork)
pca.plot + h1+plot_layout(guides = "collect")
ggsave("plot/DESeq2_heat_ve_pca.png",width = 10,height = 7)
#write.csv(dat,file = 'export/tcga/auc_second_tcga.csv')
save(dat,cg1,file = 'Rdata/tcga/pre/exp_tcga_deg.Rdata')
