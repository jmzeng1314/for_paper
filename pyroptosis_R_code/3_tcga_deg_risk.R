rm(list = ls())
load('database/TCGA-LIHC_gdc.Rdata')
load('Rdata/tcga/pre/clinicaldata_juzhenshuju.Rdata')
library(tinyarray)
library(tidyverse)
exp[1:4,1:4]
exp = trans_exp(exp,mrna_only = T)
exp[1:4,1:4]
Group=ifelse(str_sub(colnames(exp),14,15)=='01','tumor','normal')
Group=factor(Group,levels = c('normal','tumor'))
exp=exp[,Group=='tumor']

exp[1:4,1:4]
dim(exp)
colnames(exp)=str_sub(colnames(exp),1,12)
exp=exp[,colnames(exp)%in%meta$ID]
dim(exp)
exp[1:4,1:4]
exprSet=log2(cpm(exp)+1)
identical(rownames(meta),colnames(exp))
Group=meta$ri
table(Group)

#deseq2----
library(DESeq2)
colData <- data.frame(row.names =colnames(exp), 
                      condition=Group)
if(!file.exists("Rdata/tcga/tcga_deg_risk_dd.Rdata")){
  dds <- DESeqDataSetFromMatrix(
    countData = exp,
    colData = colData,
    design = ~ condition)
  dds <- DESeq(dds)
  save(dds,file = "Rdata/tcga/tcga_deg_risk_dd.Rdata")
}
load(file ='Rdata/tcga/tcga_deg_risk_dd.Rdata')
class(dds)
res <- results(dds, contrast = c("condition",rev(levels(Group))))
c("condition",rev(levels(Group)))
class(res)
DEG1 <- as.data.frame(res)
DEG1 <- DEG1[order(DEG1$pvalue),] 
DEG1 = na.omit(DEG1)
head(DEG1)

#添加change列标记基因上调下调
logFC_t = 1
pvalue_t = 0.05

k1 = (DEG1$pvalue < pvalue_t)&(DEG1$log2FoldChange < -logFC_t);table(k1)

k2 = (DEG1$pvalue < pvalue_t)&(DEG1$log2FoldChange > logFC_t);table(k2)
DEG1$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG1$change)
head(DEG1)
library(ggplot2)
library(tinyarray)
exp[1:4,1:4]
# cpm 去除文库大小的影响
dat=log2(cpm(exp)+1)
a=rownames(DEG1)[DEG1$change!='NOT']
pca.plot = draw_pca(dat[a,],Group,addEllipses = F);pca.plot
save(pca.plot,file = "plot/tcga/pre/deg_risk_pcaplot.Rdata")
ggsave('plot/tcga/pre/tcga_risk_deg_pcaplot.png',height = 5,width = 5)
cg1 = rownames(DEG1)[DEG1$change !="NOT"]
identical(colnames(exprSet),colnames(exp))
identical(rownames(meta),colnames(exprSet))
table(Group)
exp_h=exprSet[a,Group=='highrisk']
exp_l=exprSet[a,Group=='lowrisk']
identical(rownames(exp_h),rownames(exp_l))
exp_dat=cbind(exp_l,exp_h)
exp_dat[1:4,1:4]
dim(exp_dat)
Group_dat=c(rep('lowrisk',228),rep('highrisk',83))
Group_dat=factor(Group_dat,levels = c('lowrisk','highrisk'))
h1 = draw_heatmap(exp_dat,Group_dat,n_cutoff = 2,cluster_cols = F)
h2=draw_heatmap(exp_dat,Group_dat,n_cutoff = 2,cluster_cols = F,annotation_legend = T);h2
ggsave('plot/tcga/pre/tcga_risk_deg_heatmap.png',height = 5,width = 5)

v1 = draw_volcano(DEG1,pkg = 1,logFC_cutoff = logFC_t);v1
ggsave('plot/tcga/pre/tcga_risk_deg_volcan.png',height = 5,width = 5)
library(patchwork)
pca.plot + h1+v1+ plot_layout(guides = "collect")
ggsave('plot/tcga/pre/tcga_deg_risk_heat_ve_pca.png',width = 15,height = 5)
save(DEG1,file = 'Rdata/tcga/pre/tcga_deg_risk.Rdata')
