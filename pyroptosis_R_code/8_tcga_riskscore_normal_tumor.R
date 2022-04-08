rm(list = ls())
load('Rdata/tcga/pre/model.Rdata')
load('database/tcga_pair_exp_clinical.Rdata')
exp_small[1:4,1:4]
exp=exp_small
exp = trans_exp(exp,mrna_only = T)
exp[1:4,1:4]
Group=ifelse(str_sub(colnames(exp),14,15)=='01','tumor','non-tumor')
Group=factor(Group,levels = c('non-tumor','tumor'))
table(Group)
exp=log2(cpm(exp)+1)
dim(exp)
exp[1:4,1:4]
exp=exp[gs,]
dat=as.data.frame(t(exp))
dat$group=Group
fp <- predict(model,newdata = dat)
dat$fp=fp
dat$riskscore=fp
names(fp)=rownames(dat)
identical(rownames(dat),names(fp))
save(dat,file = 'Rdata/tcga/pre/tumor_normal_riskscore_tcga_baohanriskscore.Rdata')
b1 = ggplot(dat = dat,aes(group,riskscore,color=group))+
  geom_boxplot()+
  geom_jitter(aes(color=group))+
  #scale_fill_manual(values = c("#1d4a9d","#e5171b"))+
  stat_compare_means()+
  theme_classic()+labs(x=' ')
ggsave('plot/tcga/pre/tcga_normal_tumor_riskscore_boxplot.png',height = 5,width = 5)
write.csv(dat,file = 'export/tcga/pre/tcga_normal_tumor_riskscore.csv')
table(dat$group)
a=dat[dat$group=='non-tumor',]
b=dat[dat$group=='tumor',]
summary(a$riskscore)
summary(b$riskscore)