rm(list = ls())
load('Rdata/tcga/pre/model.Rdata')
load('database/icgc_tumor_normal_rawcount.Rdata')
library(stringr)
Group=ifelse(str_detect(colnames(exp),'cancer'),'tumor','normal')
exp=exp[,Group=='tumor']
exp=log2(cpm(exp)+1)
exp=exp[rownames(exp)%in%gs,]

dat=exp
e=t(dat)
#colnames(e)= str_replace_all(colnames(e),"-","_")
identical(str_split(rownames(e),'-',simplify = T)[,1],rownames(cl_tumor))
meta=cl_tumor
dat=cbind(meta,e)
dat$gender=as.numeric(factor(dat$gender))
dat$stage=as.numeric(factor(dat$stage))
summary(meta$time)
colnames(dat)
fp <- predict(model,newdata = dat)
meta$fp=fp
meta$riskscore=fp
names(fp)=rownames(dat)
identical(rownames(dat),names(fp))
library(Hmisc)
options(scipen=200)
with(dat,rcorr.cens(fp,Surv(time, event)))
###绘制timeROC
colnames(meta)
dat_time = meta[,c(3,5,6)]
head(dat_time)                
library(survminer)
library(survival)
library(timeROC)
result <-with(dat_time, timeROC(T=time,
                                delta=event,
                                marker=fp,
                                cause=1,
                                times=c(12,24,36),
                                iid = TRUE))
plot_dat = data.frame(fpr = as.numeric(result$FP),
                      tpr = as.numeric(result$TP),
                      time = rep(as.factor(c(12,24,36)),each = nrow(result$TP)))
library(ggplot2)
ggplot() + 
  geom_line(data = plot_dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,2,3),"-y survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
ggsave('plot/icgc/pre/cox__timeROC_icgc.png',height = 5,width = 5)
###绘制生存曲线
res.cut = survminer::surv_cutpoint(meta, time = "time", event = "event", variables = "riskscore")
cut = res.cut[["cutpoint"]][1, 1]
cut
png('plot/icgc/pre/cut_tcga.png',height = 500,width = 500)
plot(res.cut, "riskscore", palette = "npg")
dev.off()

ri = ifelse(fp<cut,"lowrisk","highrisk")
ri = factor(ri,levels = c("lowrisk","highrisk"))
names(ri)=names(fp)
meta$ri=ri
table(meta$ri)
sfit <- survfit(Surv(time, event)~ri, data=meta)
summary(sfit)$table
p_sur_train=ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
                       risk.table =TRUE,pval =TRUE,
                       conf.int =TRUE,xlab ="Time in months", 
                       ggtheme =theme_light(), 
                       ncensor.plot = TRUE)
library(patchwork)
p1=p_sur_train$plot
p2=p_sur_train$table
p3=p_sur_train$ncensor.plot
p1 /p2  + plot_layout(design = 'A
                         A
                         A
                         A
                         B')

ggsave(file='plot/icgc/pre/icgc_cox_model_kmplot.png',height = 5,width = 5)

##绘制箱线图
library(tinyarray)
identical(str_split(colnames(exp),'-',simplify = T)[,1],rownames(meta))
draw_boxplot(exp[gs,],ri,color = c("#1d4a9d","#e5171b"))
ggsave(file='plot/icgc/pre/cox_6gene_boxplot.png',height = 5,width = 5)
###绘制三图联动
fp_dat=data.frame(patientid=1:length(fp),
                  fp=as.numeric(sort(fp)),
                  ri= ri[order(fp)])

sur_dat=data.frame(patientid=1:length(fp),
                   time=meta[order(fp),'time'] ,
                   event=meta[order(fp),'event']) 


sur_dat$event=ifelse(sur_dat$event==0,'alive','death')
sur_dat$event=factor(sur_dat$event,levels = c("death","alive"))
library(stringr)


table(gs%in%colnames(e))
a=e[order(fp),]
b=a[,colnames(a)%in%gs]
colnames(b)
exp_dat=t(b)

colnames(exp_dat) = str_sub(colnames(exp_dat),1,12)
rownames(exp_dat) = str_replace_all(rownames(exp_dat),"_","-")
###第一个图----
p1=ggplot(fp_dat,aes(x=patientid,y=fp))+
  geom_point(aes(color = ri))+
  scale_color_manual(values = c("blue","red"))+
  geom_vline(xintercept = sum(fp<cut),lty = 2)+
  scale_x_continuous(expand=c(0.01,0))+
  labs(x = "",y = "risk score")+
  theme_bw()
#第二个图
p2=ggplot(sur_dat,aes(x=patientid,y=time))+
  geom_point(aes(col=event))+
  scale_color_manual(values = c("red","blue"))+
  geom_vline(xintercept = sum(fp<cut),lty = 2)+
  scale_x_continuous(expand=c(0.01,0))+
  labs(x = "")+
  theme_bw()
#第三个图
a=c('TP53','GPX4','GSDMC','BAK1')
exp_dat=exp_dat[match(a,rownames(exp_dat)),]
annotation_col = data.frame(group = ri,
                            row.names = names(ri))
mycolors <- colorRampPalette(c("blue","white", "red"), bias = 1.2)(100)
ann_colors = list(
  group = c(lowrisk="blue", highrisk="red")
)
p3=pheatmap::pheatmap(exp_dat,
                      col= mycolors,
                      annotation_col = annotation_col,
                      annotation_colors =ann_colors,
                      scale = "row",
                      breaks = seq(-3,3,length.out = 100),
                      show_colnames = F,
                      cluster_cols = F,
                      annotation_legend = F)

library(tinyarray)
n = scale(t(exp_dat))
n[n > 3] = 3
n[n < -3] = -3
p3 = ggheat(n,fp_dat$ri,show_rownames = F,legend_color= c("blue","red"),color = c("blue","white","red"))+
  theme(axis.text = element_text(size = 8))
p3
colnames(n)
#拼图实现三图联动
library(patchwork)
p1 /p2 /p3 + plot_layout(design = 'A
                         B
                         B
                         C
                         C
                         C')
ggsave("plot/icgc/pre/fengxianyinzi_santuliandong_lass_coxzhubuduoyuanhuigui.png",height = 5,width = 5)
save(meta,exp_dat,file = 'Rdata/icgc/pre/linchuangzuiquan_icgc.Rdata')

rm(list = ls())
load('Rdata/icgc/pre/linchuangzuiquan_icgc.Rdata')
meta=meta[match(str_split(colnames(exp_dat),'-',simplify = T)[,1],rownames(meta)),]
identical(str_split(colnames(exp_dat),'-',simplify = T)[,1],rownames(meta))
#ri=ri[match(rownames(meta),names(ri))]
#identical(names(ri),rownames(meta))
#identical(names(ri),str_split(colnames(exp_dat),'-',simplify = T)[,1])
meta$age_group=ifelse(meta$age<65,'younger','older')
meta$age_group=factor(meta$age_group,levels = c('younger','older'))
table(meta$age_group)

# colnames(dat)
# identical(rownames(dat),rownames(meta))
# dat=dat[match(rownames(meta),rownames(dat)),]
# identical(rownames(dat),rownames(meta))
# dat$ri=meta$ri
# colnames(dat)
# colnames(meta)

meta_yonuger=meta[meta$age_group=='younger',]
sfit <- survfit(Surv(time, event)~ri, data=meta_yonuger)
summary(sfit)$table
p_younger=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),title='Younger')

meta_older=meta[meta$age_group=='older',]
sfit <- survfit(Surv(time, event)~ri, data=meta_older)
summary(sfit)$table
p_older=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),title='Older')

table(meta$gender)
meta_femal=meta[meta$gender=='female',]
sfit <- survfit(Surv(time, event)~ri, data=meta_femal)
summary(sfit)$table
p_femal=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),title='Female')

meta_mal=meta[meta$gender=='male',]
sfit <- survfit(Surv(time, event)~ri, data=meta_mal)
summary(sfit)$table
p_mal=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Male')



meta_early_stage=meta[meta$stage%in% c(1,2),]
sfit <- survfit(Surv(time, event)~ri, data=meta_early_stage)
summary(sfit)$table
p_early_stage=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Stage: I-II')

meta_later_stage=meta[meta$stage%in% c(3,4),]
sfit <- survfit(Surv(time, event)~ri, data=meta_later_stage)
summary(sfit)$table
p_later_stage=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Stage: III-IV')

save(p_younger,p_older,p_mal,p_femal,p_early_stage,p_later_stage,file = 'Rdata/icgc/pre/icgc_km_clinical.Rdata')
dev.off()
library(patchwork)
(p_younger$plot|p_older$plot|p_mal$plot)/
  (p_femal$plot|p_early_stage$plot|p_later_stage$plot)
ggsave('plot/icgc/butongzu_risk_kmplot.png',width = 20,height = 10)
str(meta)
meta=meta[!is.na(meta$age),]
table(meta$age_group)
a=meta[meta$age_group=='younger',]
b=meta[meta$age_group=='older',]
summary(a$riskscore)
summary(b$riskscore)
b_age = ggplot(dat = meta,aes(age_group,fp))+
  geom_boxplot(aes(fill = age_group))+
  scale_fill_manual(values = c("#1d4a9d","#e5171b"))+
  stat_compare_means(label.y = 3)+
  theme_classic()+labs(x=' ',y='risk score')
table(meta$gender)
a=meta[meta$gender=='male',]
b=meta[meta$gender=='female',]
summary(a$riskscore)
summary(b$riskscore)
b_gender=ggplot(dat = meta,aes(gender,fp))+
  geom_boxplot(aes(fill = gender))+
  stat_compare_means(label.y = 3)+
  scale_fill_manual(values = c("#1d4a9d","#e5171b"))+
  stat_compare_means(method = "wilcox.test",comparisons =list(c('MALE','FEMALE')))+
  theme_classic()+labs(x=' ',y='risk score')

table(meta$stage)
meta$stage=ifelse(meta$stage==1,'I',
                  ifelse(meta$stage==2,'II',
                         ifelse(meta$stage==3,'III','IV')))
meta$stage=as.factor(meta$stage)
a=meta[meta$stage=='I',]
b=meta[meta$stage=='II',]
e=meta[meta$stage=='III',]
f=meta[meta$stage=='IV',]
summary(a$riskscore)
summary(b$riskscore)
summary(e$riskscore)
summary(f$riskscore)

b_stage=ggplot(dat = meta,aes(stage,fp))+
  geom_boxplot(aes(fill = stage))+
  scale_fill_manual(values = c("#1d4a9d","#9F5F9F","#DB7093","#e5171b"))+
  stat_compare_means(label.y = 5)+
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("I","II"),
                                        c("I","III"),
                                        c("I","IV"),
                                        c('II','III'),
                                        c('II','IV'),
                                        c('III','IV')),
                     #symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
  )+
  theme_classic()+labs(x=' ',y='risk score')
save(b_age,b_gender,b_stage,file='Rdata/icgc/pre/icgc_boxplot_clinical.Rdata')
library(patchwork)
(b_age|b_gender|b_stage)
ggsave('plot/icgc/boxplot_butongzu.png',height = 4,width = 15)



##1.age
m = coxph(Surv(time, event) ~ age, data =  meta)
beta <- coef(m)
se <- sqrt(diag(vcov(m)))
HR <- exp(beta)
HRse <- HR * se
tmp_age <- round(cbind(coef = beta,
                       se = se, z = beta/se,
                       p = 1 - pchisq((beta/se)^2, 1),
                       HR = HR, HRse = HRse,
                       HRz = (HR - 1) / HRse,
                       HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                       HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                       HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
###2.gender
m = coxph(Surv(time, event) ~ gender, data =  meta)
beta <- coef(m)
se <- sqrt(diag(vcov(m)))
HR <- exp(beta)
HRse <- HR * se
tmp_gender <- round(cbind(coef = beta,
                          se = se, z = beta/se,
                          p = 1 - pchisq((beta/se)^2, 1),
                          HR = HR, HRse = HRse,
                          HRz = (HR - 1) / HRse,
                          HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                          HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                          HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)

colnames(dat)
meta$stage=as.numeric(factor(meta$stage))
table(meta$stage)
##3.stage
m = coxph(Surv(time, event) ~stage , data =  meta)
beta <- coef(m)
se <- sqrt(diag(vcov(m)))
HR <- exp(beta)
HRse <- HR * se
tmp_stage <- round(cbind(coef = beta,
                         se = se, z = beta/se,
                         p = 1 - pchisq((beta/se)^2, 1),
                         HR = HR, HRse = HRse,
                         HRz = (HR - 1) / HRse,
                         HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                         HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                         HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
colnames(dat)

##9.riskscore

m = coxph(Surv(time, event) ~riskscore , data =  meta)
beta <- coef(m)
se <- sqrt(diag(vcov(m)))
HR <- exp(beta)
HRse <- HR * se
tmp_riskscore <- round(cbind(coef = beta,
                             se = se, z = beta/se,
                             p = 1 - pchisq((beta/se)^2, 1),
                             HR = HR, HRse = HRse,
                             HRz = (HR - 1) / HRse,
                             HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                             HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                             HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)

tmp=rbind(tmp_age,tmp_gender)
tmp=rbind(tmp,tmp_stage)
tmp=rbind(tmp,tmp_riskscore)

save(tmp,file='Rdata/icgc/uni_cox_clinical.Rdata')

write.csv(tmp,file = 'export/icgc_uni_cox_clinical.csv')
rownames(tmp)
s = as.formula(paste("Surv(time, event) ~ ",paste(c(
  "age",
  "gender",
  "stage",
  "riskscore"),collapse = "+")))
s
model <- coxph(s, data = meta )
summary(model)$concordance
dev.off()
png('plot/icgc/all_clinical_forest.png',height = 400,width = 400)
ggforest(model, data =meta)
dev.off()
save(model,s,file = 'Rdata/icgc/clinical_cox.Rdata')

rm(list=ls())
load('Rdata/icgc/clinical_cox.Rdata')
##多因素cox风险森林图

m = summary(model)
colnames(m$coefficients)
colnames(m$conf.int)
#p值改一下格式，加上显著性
p = ifelse(
  m$coefficients[, 5] < 0.001,
  "<0.001 ***",
  ifelse(
    m$coefficients[, 5] < 0.01,
    "<0.01  **",
    ifelse(
      m$coefficients[, 5] < 0.05,
      paste(round(m$coefficients[, 5], 3), " *"),
      round(m$coefficients[, 5], 3)
    )
  )
)
#HR和它的置信区间
dat2 = as.data.frame(round(m$conf.int[, c(1, 3, 4)], 2))
dat2 = tibble::rownames_to_column(dat2, var = "Trait")
colnames(dat2)[2:4] = c("HR", "lower", "upper")
dat1=m$coefficients
dat2$Coef=dat1[,1]
#需要在图上显示的HR文字和p值
dat2$HR2 = paste0(dat2[, 2], "(", dat2[, 3], "-", dat2[, 4], ")")
dat2$p = p

str(dat2)

dat2$Trait=c('age',"gender",
             "stage", 
             'riskscore')

dat1=m$coefficients

library(survival)
library(survminer)
library(forestplot)
library(stringr)
dat2 = rbind(
  c("Trait", NA, NA, NA,"Coef", "HR", "P"),
  dat2[1, ],
  dat2[2:nrow(dat2), ]
)
dat2$Coef=c("Coef","0.001","0.966", "0.787","0.784")
str(dat2)
dat2$HR=as.numeric(dat2$HR)
dat2$lower=as.numeric(dat2$lower)
dat2$upper=as.numeric(dat2$upper)

colnames(dat2)
png('plot/icgc/multi_clinical_cox.png',height = 300,width =500 )
forestplot(
  dat2[, c(1, 5, 6,7)],
  mean = dat2[, 2],
  lower = dat2[, 3],
  upper = dat2[, 4],
  zero = 1,
  boxsize = 0.2,
  col = fpColors(box = '#1075BB', lines = 'black', zero = 'grey'),
  lty.ci = "solid",
  graph.pos = 3,
  #xticks = F,
  is.summary = c(T, rep(F, 4)),
  align = "l",
  hrzl_lines = list(
    "1" = gpar(lty=1.5),
    "2" = gpar(lty=1.5),
    "6"= gpar(lty=1.5)),
  colgap = unit(5, 'mm')
)
dev.off()


###单因素
rm(list=ls())
uni_clinical=read.csv('export/icgc/pre/icgc_uni_cox_clinical_zhengli.csv')
tmp=uni_clinical
tmp$X=c('age','gender','stage','riskscore')
tmp=rbind(c("Trait","Coef" ,NA, NA, NA, "HR", "P"),tmp)
str(tmp)
tmp$HR=as.numeric(tmp$HR)
tmp$lower=as.numeric(tmp$lower)
tmp$upper=as.numeric(tmp$upper)
dev.off()
library(forestplot)
png('plot/icgc/pre/unni_clinical_cox.png',height = 300,width =500 )
forestplot(
  tmp[, c(1, 2, 6,7)],
  mean = tmp[, 3],
  lower = tmp[, 4],
  upper = tmp[, 5],
  zero = 1,
  boxsize = 0.2,
  col = fpColors(box = '#1075BB', lines = 'black', zero = 'grey'),
  lty.ci = "solid",
  graph.pos = 3,
  #xticks = F,
  is.summary = c(T, rep(F, 4)),
  align = "l",
  hrzl_lines = list(
    "1" = gpar(lty=1.5),
    "2" = gpar(lty=1.5),
    "6"= gpar(lty=1.5)),
  colgap = unit(5, 'mm')
)
dev.off()

rm(list = ls())
load('Rdata/icgc/pre/linchuangzuiquan_icgc.Rdata')
colnames(exp_dat)=str_split(colnames(exp_dat),'-',simplify = T)[,1]
identical(rownames(meta),colnames(exp_dat))
exp_dat=exp_dat[,match(rownames(meta),colnames(exp_dat))]
identical(rownames(meta),colnames(exp_dat))
Group=meta$ri
table(Group)
Group=factor(Group,levels = c('lowrisk','highrisk'))
pca.plot = draw_pca(exp_dat,Group,addEllipses = F);pca.plot
ggsave('plot/icgc/pre/pca_risk.png',height = 5,width = 5)
