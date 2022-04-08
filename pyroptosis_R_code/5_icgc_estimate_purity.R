rm(list = ls())
load('database/icgc_tumor_normal_rawcount.Rdata')
library(stringr)
Group=ifelse(str_detect(colnames(exp),'cancer'),'tumor','normal')
exp=exp[,Group=='tumor']
colnames(exp)=str_split(colnames(exp),'-',simplify = T)[,1]
load('Rdata/icgc/pre/linchuangzuiquan_icgc.Rdata')
identical(rownames(meta),colnames(exp))
dim(exp)
colnames(exp)
exp[1:4,1:4]
if(!require("estimate")){
  library(utils)
  rforge <- "http://r-forge.r-project.org"
  install.packages("estimate", repos=rforge, dependencies=TRUE)
}
library(estimate)
library(stringr)
exp[1:4,1:4]
dat = log2(edgeR::cpm(exp)+1)

estimate <- function(dat,pro){
  input.f='Rdata/icgc/pre/icgc_estimate_input.txt'
  output.f='Rdata/icgc/pre/icgc_estimate_gene.gct'
  output.ds='Rdata/icgc/pre/icgc_estimate_score.gct'
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")   ## 注意platform
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
scores=estimate(dat)
scores = rownames_to_column(as.data.frame(scores),"id")
scores$id = str_replace_all(scores$id,"\\.","-")
head(scores)
TumorPurity = cos(0.6049872018+0.0001467884 * scores[,3])
head(TumorPurity)
scores$TumorPurity=TumorPurity
dim(scores)
identical(scores$id,rownames(meta))
a=scores
class(a)
a=as.data.frame(a)
a$Group=meta$ri
#scores$Group =cl$ri
library(ggplot2)
library(ggpubr)
b1 = ggplot(dat = a,aes(Group,StromalScore))+
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#1d4a9d","#e5171b"))+
  stat_compare_means()+
  theme_classic()+labs(x=' ')
b2 = ggplot(dat = a,aes(Group,ImmuneScore ))+
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#1d4a9d","#e5171b"))+
  stat_compare_means()+
  theme_classic()+labs(x=' ')
colnames(a)
b3=ggplot(dat = a,aes(Group, ESTIMATEScore))+
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#1d4a9d","#e5171b"))+
  stat_compare_means()+
  theme_classic()+labs(x=' ')
b4=ggplot(dat = a,aes(Group, TumorPurity))+
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#1d4a9d","#e5171b"))+
  stat_compare_means()+
  theme_classic()+labs(x=' ')
b1|b2|b3|b4
ggsave('plot/icgc/pre/icgc_estimate.png',height = 5,width = 16)
dev.off()
b1
ggsave('plot/icgc/pre/icgc_estimate_StromalScore.png',height = 5,width = 5)
dev.off()
b2
ggsave('plot/icgc/pre/icgc_estimate_ImmuneScore.png',height = 5,width = 5)
dev.off()
b3
ggsave('plot/icgc/pre/icgc_estimate_ESTIMATEScore.png',height = 5,width = 5)
dev.off()
b4
ggsave('plot/icgc/pre/icgc_estimate_TumorPurity.png',height = 5,width = 5)
dev.off()
table(a$Group)
b=a[a$Group=='lowrisk',]
e=a[a$Group=='highrisk',]
summary(b$StromalScore)
summary(b$ImmuneScore)
summary(b$ESTIMATEScore)
summary(b$TumorPurity)
summary(e$StromalScore)
summary(e$ImmuneScore)
summary(e$ESTIMATEScore)
summary(e$TumorPurity)
