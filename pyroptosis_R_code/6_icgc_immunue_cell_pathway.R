rm(list = ls())
library(GSVA)
##绘制immunecell
load('database/icgc_tumor_normal_rawcount.Rdata')
library(stringr)
Group=ifelse(str_detect(colnames(exp),'cancer'),'tumor','normal')
exp=exp[,Group=='tumor']
colnames(exp)=str_split(colnames(exp),'-',simplify = T)[,1]

dim(exp)
exp[1:4,1:4]
identical(colnames(exp),rownames(cl_tumor))
geneset=readxl::read_xlsx('database/immune_cell_16.xlsx')
geneset = split(geneset$Markers,geneset$`Immune cells/pathways`)
lapply(geneset[1:3], head)
class(exp)
exp=as.matrix(exp)
f = "Rdata/icgc/pre/icgc_risk_deg_immu_cell16.Rdata"
if(!file.exists(f)){
  re <- gsva(exp, geneset, method="ssgsea",
             mx.diff=FALSE, verbose=FALSE,kcdf=c('Poisson'))
  save(re,file = f)
}
load('Rdata/icgc/pre/linchuangzuiquan_icgc.Rdata')
identical(rownames(meta),colnames(exp))
load(f)
cl=meta
table(cl$ri)
colnames(cl)
class(cl$ri)
identical(colnames(exp),rownames(cl))
Group=cl$ri
draw_boxplot(re,Group,color = c("#1d4a9d","#e5171b"),xlab = ' ',ylab = 'Score',sort = F)
ggsave('plot/icgc/pre/icgc_ssgsva_immune_cell_13_deg_risk.png',height = 5,width = 8)


##绘制pathway的通路图

wayset=readxl::read_xlsx('database/immune_pathway_13.xlsx')
wayset = split(wayset$Markers,wayset$`Immune cells/pathways`)
lapply(wayset[1:3], head)

class(exp)
exp=as.matrix(exp)
identical(colnames(exp),rownames(cl))
e = "Rdata/icgc/pre/icgc_immu_pathway13_deg_risk.Rdata"
if(!file.exists(e)){
  re <- gsva(exp, wayset, method="ssgsea",
             mx.diff=FALSE, verbose=FALSE,kcdf=c('Poisson'))
  save(re,file = e)
}

load(e)
colnames(cl)
class(cl$ri)
identical(colnames(exp),rownames(cl))
Group=cl$ri
table(Group)
draw_boxplot(re,Group,color = c("#1d4a9d","#e5171b"),xlab = ' ',ylab = 'Score',sort = F)
ggsave('plot/icgc/pre/icgc_ssgsva_immune_pathway_13.png',height = 5,width = 8)




