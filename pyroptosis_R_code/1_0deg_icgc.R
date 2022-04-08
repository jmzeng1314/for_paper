rm(list = ls())
load('database/icgc_tumor_normal_rawcount.Rdata')
colnames(exp)
exp_normal=exp[,str_detect(colnames(exp),'normal')]
str_split(colnames(exp_normal),'-',simplify = T)[,1]
patient=str_split(colnames(exp_normal),'-',simplify = T)[,1]
cl_pair=cl_tumor[rownames(cl_tumor)%in%patient,]
k=str_split(colnames(exp),'-',simplify = T)[,1]
exp=exp[,k%in%patient]
dim(exp)
Group=ifelse(str_detect(colnames(exp),'normal'),'non-tumor','tumor')
Group=factor(Group,levels = c('non-tumor','tumor'))
save(exp,Group,cl_pair,file = 'export/icgc/pair_icgc_tumor_normal_rawdata.Rdata')
a=exp
colnames(exp)
colnames(Group)
load('database/pyroptosis.Rdata')
exp=exp[rownames(exp)%in%pyroptosis,]
###做差异分析
class(exp)

library(DESeq2)
colData <- data.frame(row.names =colnames(exp), 
                      condition=Group)
dim(exp)

if(!file.exists("Rdata/pyroptosis_DESeq_icgc.Rdata")){
  dds <- DESeqDataSetFromMatrix(
    countData = exp,
    colData = colData,
    design = ~ condition)
  dds <- DESeq(dds)
  save(dds,file = "Rdata/pyroptosis_DESeq_icgc.Rdata")
}
load(file ="Rdata/pyroptosis_DESeq_icgc.Rdata")
class(dds)
res <- results(dds, contrast = c("condition",rev(levels(Group))))
c("condition",rev(levels(Group)))

class(res)

DEG1 <- as.data.frame(res)
DEG1 <- DEG1[order(DEG1$pvalue),] 
DEG1 = na.omit(DEG1)
head(DEG1)


#添加change列标记基因上调下调
logFC_t =0
pvalue_t = 0.05

k1 = (DEG1$padj < pvalue_t)&(DEG1$log2FoldChange < -logFC_t);table(k1)
k2 = (DEG1$padj < pvalue_t)&(DEG1$log2FoldChange > logFC_t);table(k2)
DEG1$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG1$change)
head(DEG1)
save(DEG1,Group,file = "Rdata/DEG1_icgc.Rdata")


write.csv(DEG1,file = 'export/deseq2_deg_icgc.csv')

#可视化
library(ggplot2)
library(tinyarray)
exp=a
exp[1:4,1:4]
# cpm 去除文库大小的影响
dat = log2(cpm(exp)+1)
dim(dat)
deg=rownames(DEG1)[DEG1$change!='NOT']
dat=dat[rownames(dat)%in%deg,]
pca.plot = draw_pca(dat,Group,addEllipses = F);pca.plot
save(pca.plot,file = "Rdata/tumor_normal_pcaplot.Rdata")
ggsave('plot/icgc_pair_pca.png',height = 5,width = 5)
dev.off()
cg1 = rownames(DEG1)[DEG1$change !="NOT"]


h1 = draw_heatmap(dat[cg1,],Group,n_cutoff = 2,cluster_cols = F,show_rownames = T)
h2 = draw_heatmap(dat[cg1,],Group,n_cutoff = 2,cluster_cols = F,show_rownames = T,annotation_legend = T,legend = T);h2
ggsave('plot/icgc_pair_heatmap.png',height = 5,width = 5)

v1 = draw_volcano(DEG1,pkg = 1,logFC_cutoff = logFC_t)




#维恩图拼图，终于搞定

library(patchwork)
#up.plot + down.plot
# 拼图
pca.plot + h1+ plot_layout(guides = "collect")
ggsave('plot/icgc_pair_heat_pca.png',width = 10,height = 7)

write.csv(dat,file = 'export/auc_second_icgc.csv')
save(dat,cg1,file ='Rdata/icgc/auc_send_icgc.Rdata')








