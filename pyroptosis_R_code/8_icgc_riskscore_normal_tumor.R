rm(list = ls())
load('Rdata/tcga/pre/model.Rdata')
load('database/icgc_pair_exp_clinical.Rdata')
exp_small[1:4,1:4]
exp=exp_small
Group=ifelse(str_split(colnames(exp),'-',simplify = T)[,2]=='cancer','tumor','non-tumor')
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
save(dat,file = 'Rdata/icgc/pre/tumor_normal_riskscore_icgc_baohanriskscore.Rdata')
b1 = ggplot(dat = dat,aes(group,riskscore,color=group))+
  geom_boxplot()+
  geom_jitter(aes(color=group))+
  #scale_fill_manual(values = c("#1d4a9d","#e5171b"))+
  stat_compare_means()+
  theme_classic()+labs(x=' ')
ggsave('plot/icgc/pre/icgc_normal_tumor_riskscore_boxplot.png',height = 5,width = 5)
write.csv(dat,file = 'export/icgc/pre/icgc_normal_tumor_riskscore.csv')
a=dat[dat$group=='non-tumor',]
b=dat[dat$group=='tumor',]
summary(a$riskscore)
summary(b$riskscore)