rm(list = ls())
library(GSVA)
##绘制immunecell
load('Rdata/tcga/pre/clinicaldata_juzhenshuju.Rdata')
meta_all=meta
load("database/TCGA-LIHC_gdc.Rdata")
library(tinyarray)
library(tidyverse)
exp[1:4,1:4]
exp = trans_exp(exp,mrna_only = T)
exp[1:4,1:4]
exp=exp[,Group=='tumor']
colnames(exp)=str_sub(colnames(exp),1,12)
exp=exp[,colnames(exp)%in%rownames(meta)]
dim(exp)
exp[1:4,1:4]
identical(colnames(exp),rownames(meta))
#exp=log2(cpm(exp)+1)
exp[1:4,1:4]
geneset=readxl::read_xlsx('database/immune_cell_16.xlsx')
geneset = split(geneset$Markers,geneset$`Immune cells/pathways`)
lapply(geneset[1:3], head)
class(exp)
exp=as.matrix(exp)####这个一定要是矩阵，不能是数据框
f = "Rdata/tcga/pre/tcga_immu_cell16.Rdata"
if(!file.exists(f)){
  re <- gsva(exp, geneset, method="ssgsea",
             mx.diff=FALSE, verbose=FALSE,,kcdf=c('Poisson'))
  save(re,file = f)
}

load(f)
table(meta$ri)
colnames(meta)
class(meta$ri)
identical(colnames(exp),rownames(meta))
Group=meta$ri
draw_boxplot(re,Group,color = c("#1d4a9d","#e5171b"),xlab = ' ',ylab = 'Score',sort = F)
ggsave('plot/tcga/pre/deg_risk_ssgsva_immune_cell_16.png',height = 5,width = 8)

##绘制pathway的通路图

wayset=readxl::read_xlsx('database/immune_pathway_13.xlsx')
wayset = split(wayset$Markers,wayset$`Immune cells/pathways`)
lapply(wayset[1:3], head)

class(exp)
exp=as.matrix(exp)
identical(colnames(exp),rownames(meta))
e = "Rdata/tcga/pre/tcga_immu_pathway13.Rdata"
if(!file.exists(e)){
  re <- gsva(exp, wayset, method="ssgsea",
             mx.diff=FALSE, verbose=FALSE,kcdf=c('Poisson'))
  save(re,file = e)
}

load(e)
colnames(meta)
class(meta$ri)
identical(colnames(exp),rownames(meta))
Group=meta$ri
draw_boxplot(re,Group,color = c("#1d4a9d","#e5171b"),xlab = ' ',ylab = 'Score',sort = F)
ggsave('plot/tcga/pre/tcga_deg_risk_ssgsva_immune_pathway_13.png',height = 5,width = 8)





