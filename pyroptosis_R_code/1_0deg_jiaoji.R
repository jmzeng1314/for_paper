rm(list = ls())
load('database/tcga_pair_exp_clinical.Rdata')
dim(exp_small)
exp=exp_small
dat=log2(cpm(exp)+1)
dat=trans_exp(dat,mrna_only = T)
dat_tcga=dat
load('Rdata/tcga/pre/exp_tcga_deg.Rdata')
dat_tcga2=dat_tcga[rownames(dat_tcga)%in%cg1,]
deg_tcga=cg1

load('database/icgc_pair_exp_clinical.Rdata')
exp=exp_small
dim(exp)
exp[1:4,1:4]
dat=log2(cpm(exp)+1)
dat_icgc=dat
load('Rdata/icgc/auc_send_icgc.Rdata')
dat_icgc=dat_icgc[rownames(dat_icgc)%in%cg1,]
deg_icgc=cg1
write.table(deg_tcga,file = 'export/tcga/deg_tcga_genename.txt',quote = F,row.names = F)
write.table(deg_icgc,file = 'export/icgc/deg_icgc_genename.txt',quote = F,row.names = F)
load('Rdata/tcga/DEG1.Rdata')
DEG1_tcga=DEG1
load('Rdata/DEG1_icgc.Rdata')
DEG1_icgc=DEG1
UP=function(df){
  rownames(df)[df$change=="UP"]
}
DOWN=function(df){
  rownames(df)[df$change=="DOWN"]
}
up=intersect(UP(DEG1_tcga),UP(DEG1_icgc))
down=intersect(DOWN(DEG1_tcga),DOWN(DEG1_icgc))
deg=c(up,down)
dat_tcga_2=dat_tcga[rownames(dat_tcga)%in%deg,]
dat_icgc_2=dat_icgc[rownames(dat_icgc)%in%deg,]
dat_tcga_2=as.data.frame(t(dat_tcga_2))
dat_tcga_2$group=ifelse(str_sub(rownames(dat_tcga_2),14,15)=="01",1,0)
table(dat_tcga_2$group)
dat_icgc_2=as.data.frame(t(dat_icgc_2))
dat_icgc_2$group=ifelse(str_detect(rownames(dat_icgc_2),'cancer'),1,0)
table(dat_icgc_2$group)

write.csv(dat_icgc_2,file = 'export/icgc/pre/deg_tumor_normal_roc_icgc.csv')
write.csv(dat_tcga_2,file = 'export/tcga/pre/deg_tumor_normal_roc_tcga.csv')
