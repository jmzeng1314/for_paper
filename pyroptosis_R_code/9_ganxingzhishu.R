rm(list=ls())
mR=read.csv('database/mRNAsi_tcga.csv',row.names = 1)
class(mR$cancer.type)
mR1=mR[mR$cancer.type%in%'LIHC',]
a=rownames(mR1)[str_sub(rownames(mR1),14,15)=='01']
mR1=mR1[rownames(mR1)%in%a,]
load('Rdata/tcga/pre/clinicaldata_juzhenshuju.Rdata')
mR1=mR1[str_sub(rownames(mR1),1,12)%in%rownames(meta),]
rownames(mR1)=str_sub(rownames(mR1),1,12)
meta=meta[rownames(meta)%in%rownames(mR1),]
identical(rownames(meta),rownames(mR1))
mR1=mR1[match(rownames(meta),rownames(mR1)),]
identical(rownames(meta),rownames(mR1))
meta$mRNAsi=mR1$mRNAsi
#检测是否符合正态分布
library(lattice)
library(MASS)
histogram(meta$mRNAsi)
histogram(meta$riskscore)#不符合正态分布
###所以用spearman相关分析进行分析
library(ggplot2)
library(ggpubr)
library(ggpubr)
ggplot(data=meta, aes(x=riskscore, y=mRNAsi))+
  geom_point(color="red")+
  stat_smooth(method="lm",se=FALSE)+stat_cor(data=meta, method = 'spearman')+
  theme_bw()
#stat_cor(data=dat, method = "pearson")意为用pearson相关进行相关性分析，可以自行更改方法
ggsave('plot/tcga/pre/tcga_mRNAsi_riskscore.png',height = 5,width = 5)
meta$ri
b1 = ggplot(dat = meta,aes(ri,mRNAsi,color=ri))+
  geom_boxplot()+
  geom_jitter(aes(color=ri))+
  #scale_fill_manual(values = c("#1d4a9d","#e5171b"))+
  stat_compare_means()+
  theme_classic()+labs(x=' ')
ggsave('plot/tcga/pre/tcga_mRNAsi_rigroup.png',height = 5,width = 5)
a=meta[meta$ri=='lowrisk',]
b=meta[meta$ri=='highrisk',]
summary(a$mRNAsi)
summary(b$mRNAsi)



###mDNAsi
rm(list = ls())
mR=read.csv('database/mDNAsi_tcga.csv',row.names = 1)
class(mR$cancer.type)
mR1=mR[mR$cancer.type%in%'LIHC',]
a=rownames(mR1)[str_sub(rownames(mR1),14,15)=='01']
mR1=mR1[rownames(mR1)%in%a,]
load('Rdata/tcga/pre/clinicaldata_juzhenshuju.Rdata')
mR1=mR1[str_sub(rownames(mR1),1,12)%in%rownames(meta),]
rownames(mR1)=str_sub(rownames(mR1),1,12)
meta=meta[rownames(meta)%in%rownames(mR1),]
identical(rownames(meta),rownames(mR1))
mR1=mR1[match(rownames(meta),rownames(mR1)),]
identical(rownames(meta),rownames(mR1))
meta$mDNAsi=mR1$mDNAsi
#检测是否符合正态分布
library(lattice)
library(MASS)
histogram(meta$mDNAsi)##不符合正态分布
histogram(meta$riskscore)#不符合正态分布
###所以用spearman相关分析进行分析
library(ggplot2)
library(ggpubr)
library(ggpubr)
ggplot(data=meta, aes(x=riskscore, y=mDNAsi))+
  geom_point(color="red")+
  stat_smooth(method="lm",se=FALSE)+stat_cor(data=meta, method = 'spearman')+
  theme_bw()
#stat_cor(data=dat, method = "pearson")意为用pearson相关进行相关性分析，可以自行更改方法
ggsave('plot/tcga/pre/tcga_mDNAsi_riskscore.png',height = 5,width = 5)

b1 = ggplot(dat = meta,aes(ri,mDNAsi,color=ri))+
  geom_boxplot()+
  geom_jitter(aes(color=ri))+
  #scale_fill_manual(values = c("#1d4a9d","#e5171b"))+
  stat_compare_means()+
  theme_classic()+labs(x=' ')
ggsave('plot/tcga/pre/tcga_mDNAsi_rigroup.png',height = 5,width = 5)
a=meta[meta$ri=='lowrisk',]
b=meta[meta$ri=='highrisk',]
summary(a$mDNAsi)
summary(b$mDNAsi)
