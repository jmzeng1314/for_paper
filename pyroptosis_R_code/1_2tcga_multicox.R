rm(list = ls())
load('Rdata/tcga/pre/cox_gene_data.Rdata')
exprSet=exp_cox_gene
dim(exprSet)
e=t(exprSet)
colnames(e)= str_replace_all(colnames(e),"-","_")
dat=cbind(meta,e)
colnames(dat)
cox
library(survival)
library(survminer)
 library(My.stepwise)
colnames(dat)
vl <- colnames(dat)[c(25:ncol(dat))]
My.stepwise.coxph(Time = "time",
                  Status = "event",
                  variable.list = vl,
                  data = dat)
model=coxph(formula = Surv(time, event) ~ GSDMC + GPX4 + TP53 + BAK1, 
            data = dat, method = "efron")
summary(model)

gs = str_split("GSDMC + GPX4 + TP53 + BAK1"," \\+ ")[[1]]
gs
save(model,gs,file='Rdata/tcga/pre/model.Rdata')
summary(model$concordance)
table(gs%in%colnames(e))
#fp <- predict(model,newdata = dat)
library(Hmisc)
#options(scipen=200)

png('plot/tcga/pre/unti_multi_cox.png',height = 300,width = 500)
ggforest(model, data =dat)
dev.off()
fp <- predict(model,data = dat)
names(fp)=dat$ID
meta$fp=fp
meta$riskscore=fp
dat$riskscore=fp
identical(names(fp),rownames(meta))
colnames(meta)
res.cut = survminer::surv_cutpoint(dat, time = "time", event = "event", variables = "riskscore",minprop = 0.2)
cut = res.cut[["cutpoint"]][1, 1]
cut
dev.off()
png('plot/tcga/cut_tcga.png',height = 500,width = 500)
plot(res.cut, "riskscore", palette = "npg")
dev.off()
ri = ifelse(fp<cut,"lowrisk","highrisk")
ri = factor(ri,levels = c("lowrisk","highrisk"))
names(ri)=names(fp)
identical(rownames(meta),names(fp))
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

ggsave(file='plot/tcga/pre/cox_model_kmplot.png',height = 5,width = 5)
dat_time = cbind(meta[,c(1,2,20)],
                 fp = fp)
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
ggsave('plot/tcga/pre/lr_cox_timeROC.png',height = 5,width = 5)

library(tinyarray)
identical(str_sub(colnames(exprSet),1,12),rownames(meta))
draw_boxplot(exprSet[gs,],ri,color = c("#1d4a9d","#e5171b"))
ggsave(file='plot/tcga/pre/cox_3gene_boxplot.png',height = 5,width = 5)

fp_dat=data.frame(patientid=1:length(fp),
                  fp=as.numeric(sort(fp)),
                  ri= ri[order(fp)])
identical(dat$ID,names(ri))

sur_dat=data.frame(patientid=1:length(fp),
                   time=meta[order(fp),'time'] ,
                   event=meta[order(fp),'event']) 


sur_dat$event=ifelse(sur_dat$event==0,'alive','death')
sur_dat$event=factor(sur_dat$event,levels = c("death","alive"))
library(stringr)

save(gs,file='export/tcga/final_gene.Rdata')
dim(e)
length(fp)

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
ggsave("plot/tcga/pre/fengxianyinzi_santuliandong_lass_coxzhubuduoyuanhuigui.png",height = 5,width = 5)
save(meta,exp_dat,file = 'Rdata/tcga/pre/clinicaldata_juzhenshuju.Rdata')
low_ri=meta[meta$ri=='lowrisk',]
table(low_ri$event)
high_ri=meta[meta$ri=='highrisk',]
table(high_ri$event)


meta=meta[match(colnames(exp_dat),rownames(meta)),]
identical(colnames(exp_dat),rownames(meta))
ri=ri[match(rownames(meta),names(ri))]
identical(names(ri),rownames(meta))
identical(names(ri),colnames(exp_dat))

group_dat=data.frame(
  group=meta$ri,
  row.names = names(ri),
  age=meta[,'age_group'] ,
  gender=meta[,'gender'],
  stage=meta[,'stage'],
  T_stage=meta[,'T_stage'],
  N_stage=meta[,'N_stage'],
  M_stage=meta[,'M_stage']
  
) 
class(group_dat$M_stage)
group_dat$gender=factor(group_dat$gender,levels = c('MALE','FEMALE'))
group_dat$stage=factor(group_dat$stage,levels = c("I","II","III","IV"))
class(exp_dat)
exp_dat=as.data.frame(exp_dat)
p3=pheatmap::pheatmap(exp_dat,
                      col= mycolors,
                      annotation_col = group_dat,
                      annotation_colors =ann_colors,
                      scale = "row",
                      breaks = seq(-3,3,length.out = 100),
                      show_colnames = F,
                      cluster_cols = F,
                      annotation_legend = T,cluster_rows = F)
ggsave("plot/tcga/pre/4gene_age_gender_heatmap.png",p3,height = 8,width = 10)
class(meta$age)
meta$age_group
colnames(dat)
identical(rownames(dat),rownames(meta))
dat=dat[match(rownames(meta),rownames(dat)),]
identical(rownames(dat),rownames(meta))
dat$ri=meta$ri
colnames(dat)
colnames(meta)
table(meta$age_group)
meta_yonuger=meta[meta$age_group=='younger',]
sfit <- survfit(Surv(time, event)~ri, data=meta_yonuger)
summary(sfit)$table
p_younger=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Younger')

meta_older=meta[meta$age_group=='older',]
sfit <- survfit(Surv(time, event)~ri, data=meta_older)
summary(sfit)$table
p_older=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Older')

meta_femal=meta[meta$gender=='FEMALE',]
sfit <- survfit(Surv(time, event)~ri, data=meta_femal)
summary(sfit)$table
p_femal=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Female')



meta_mal=meta[meta$gender=='MALE',]
sfit <- survfit(Surv(time, event)~ri, data=meta_mal)
summary(sfit)$table
p_mal=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Male')

table(meta$family_history)
meta_hist=meta[meta$family_history=='YES',]
sfit <- survfit(Surv(time, event)~ri, data=meta_hist)
summary(sfit)$table
p_hist=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Family history: Yes')

meta_nonhist=meta[meta$family_history=='NO',]
sfit <- survfit(Surv(time, event)~ri, data=meta_nonhist)
summary(sfit)$table
p_nonhist=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Family history: No')

table(meta$inflammation_grade)
meta_mild_inflammation=meta[meta$inflammation_grade%in%c('G1','G2'),]
sfit <- survfit(Surv(time, event)~ri, data=meta_mild_inflammation)
summary(sfit)$table
p_mild_inflammation=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Inflammation_grade:G1-G2')


meta_sever_inflammation=meta[meta$inflammation_grade%in% c('G3','G4'),]
sfit <- survfit(Surv(time, event)~ri, data=meta_sever_inflammation)
summary(sfit)$table
p_sever_inflammation=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Inflammation_grade:G3-G4')

meta_early_stage=meta[meta$stage%in% c('I','II'),]
sfit <- survfit(Surv(time, event)~ri, data=meta_early_stage)
summary(sfit)$table
p_early_stage=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Stage: I-II')

meta_later_stage=meta[meta$stage%in% c('III','IV'),]
sfit <- survfit(Surv(time, event)~ri, data=meta_later_stage)
summary(sfit)$table
p_later_stage=ggsurvplot(sfit,pval=TRUE,palette = c("#E7B800", "#2E9FDF"),conf.int =F,title='Stage: III-IV')


dev.off()
library(patchwork)
  (p_mild_inflammation$plot|p_sever_inflammation$plot|p_early_stage$plot|p_later_stage$plot)
table(dat$stage)
ggsave('plot/tcga/pre/butongzu_risk_kmplot.png',width = 20,height = 5)
save(p_younger,p_older,p_mal,p_femal,p_hist,p_nonhist,
  p_mild_inflammation,p_sever_inflammation,p_early_stage,p_later_stage,file='Rdata/tcga/pre/train_tcga_km_clinical.Rdata')

meta$gender=factor(meta$gender,levels=c('MALE','FEMALE'))
meta$stage=factor(meta$stage,levels = c('I','II','III','IV'))
meta$inflammation_grade=factor(meta$inflammation_grade,
                               levels = c('G1','G2','G3','G4'))

meta$family_history=factor(meta$family_history,levels=c('NO','YES'))
meta$adjacent_inflammation=factor(meta$adjacent_inflammation,
                                  levels = c('None','Mild','Severe'))


colnames(meta)
dat=meta
dat=dat[!is.na(dat$age_group),]
table(dat$age_group)
a=dat[dat$age_group=='older',]
b=dat[dat$age_group=='younger',]
summary(a$riskscore)
summary(b$riskscore)
class(b$riskscore)
table(meta$gender)
a=meta[meta$gender=='MALE',]
summary(a$riskscore)
b=meta[meta$gender=='FEMALE',]
summary(b$riskscore)
dat$age_group=factor(dat$age_group,levels = c('younger','older'))
b_age = ggplot(dat = dat,aes(age_group,fp))+
  geom_boxplot(aes(fill = age_group))+
  scale_fill_manual(values = c("#1d4a9d","#e5171b"))+
  stat_compare_means(label.y = 3)+
  theme_classic()+labs(x=' ',y='risk score')
b_gender=ggplot(dat = meta,aes(gender,fp))+
  geom_boxplot(aes(fill = gender))+
  stat_compare_means(label.y = 3)+
  scale_fill_manual(values = c("#1d4a9d","#e5171b"))+ theme_classic()+labs(x=' ',y='risk score')
dat=meta
dat=dat[!is.na(dat$family_history),]
a=dat[dat$family_history=='YES',]
b=dat[dat$family_history=='NO',]
summary(a$riskscore)
summary(b$riskscore)
b_familyhistory=ggplot(dat = dat,aes(family_history,fp))+
  geom_boxplot(aes(fill = family_history))+
  stat_compare_means(label.y = 3)+
  scale_fill_manual(values = c("#1d4a9d","#e5171b"))+theme_classic()+labs(x=' ',y='risk score')
table(meta$stage)

dat=meta
dat=dat[!is.na(dat$stage),]
a=dat[dat$stage=='I',]
b=dat[dat$stage=='II',]
e=dat[dat$stage=='III',]
f=dat[dat$stage=='IV',]
summary(a$riskscore)
summary(b$riskscore)
summary(e$riskscore)
summary(f$riskscore)

b_stage=ggplot(dat = dat,aes(stage,fp))+
  geom_boxplot(aes(fill = stage))+
  scale_fill_manual(values = c("#1d4a9d","#9F5F9F","#DB7093","#e5171b"))+
  stat_compare_means(label.y = 6)+
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
colnames(meta)
table(meta$inflammation_grade)

dat=meta
dat=dat[!is.na(dat$inflammation_grade),]
a=dat[dat$inflammation_grade=='G1',]
b=dat[dat$inflammation_grade=='G2',]
e=dat[dat$inflammation_grade=='G3',]
f=dat[dat$inflammation_grade=='G4',]
summary(a$riskscore)
summary(b$riskscore)
summary(e$riskscore)
summary(f$riskscore)
b_inflammation=ggplot(dat = dat,aes(inflammation_grade,fp))+
  geom_boxplot(aes(fill = inflammation_grade))+
  scale_fill_manual(values = c("#1d4a9d","#9F5F9F","#DB7093","#e5171b"))+
  stat_compare_means(label.y = 5)+
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("G1","G2"),
                                        c("G1","G3"),
                                        c("G1","G4"),
                                        c('G2','G3'),
                                        c('G2','G4'),
                                        c('G3','G4'))
  )+
  theme_classic()+labs(x=' ',y='risk score')

library(patchwork)
  (b_inflammation|b_stage)
ggsave('plot/tcga/pre/clinical_pathology_drawplot.png',height = 5,width = 10)
save(b_age,b_gender,b_familyhistory,b_stage,b_inflammation,file='Rdata/tcga/pre/tcga_boxplot_clinical.Rdata')

rm(list = ls())
load('Rdata/tcga/pre/clinicaldata_juzhenshuju.Rdata')
identical(rownames(meta),colnames(exp_dat))
exp_dat=exp_dat[,match(rownames(meta),colnames(exp_dat))]
identical(rownames(meta),colnames(exp_dat))
Group=meta$ri
Group=factor(Group,levels = c('lowrisk','highrisk'))
pca_sur = draw_pca(exp_dat,Group,addEllipses = F)
ggsave(pca_sur,file='plot/tcga/pre/tcga_pca_sur.png',height = 5,width = 5)

rm(list = ls())
load('Rdata/tcga/pre/clinicaldata_juzhenshuju.Rdata')
str(meta)
meta$gender=factor(meta$gender,levels = c('MALE','FEMALE'))
meta$stage=as.numeric(factor(meta$stage))
meta$inflammation_grade=as.numeric(factor(meta$inflammation_grade))
meta$family_history=factor(meta$family_history,levels=c('NO','YES'))
options(digits = 3)
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

##3.stage
meta$stage=as.numeric(factor(meta$stage))
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
##4.family_history
m = coxph(Surv(time, event) ~family_history , data =  meta)
beta <- coef(m)
se <- sqrt(diag(vcov(m)))
HR <- exp(beta)
HRse <- HR * se
tmp_family_history <- round(cbind(coef = beta,
                                  se = se, z = beta/se,
                                  p = 1 - pchisq((beta/se)^2, 1),
                                  HR = HR, HRse = HRse,
                                  HRz = (HR - 1) / HRse,
                                  HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                                  HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                                  HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)


##5.inflammation_grade
meta$inflammation_grade=as.numeric(factor(meta$inflammation_grade))
m = coxph(Surv(time, event) ~inflammation_grade , data =  meta)
beta <- coef(m)
se <- sqrt(diag(vcov(m)))
HR <- exp(beta)
HRse <- HR * se
tmp_inflammation_grade <- round(cbind(coef = beta,
                                      se = se, z = beta/se,
                                      p = 1 - pchisq((beta/se)^2, 1),
                                      HR = HR, HRse = HRse,
                                      HRz = (HR - 1) / HRse,
                                      HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                                      HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                                      HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
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
tmp=rbind(tmp,tmp_family_history)
tmp=rbind(tmp,tmp_inflammation_grade)
tmp=rbind(tmp,tmp_stage)
tmp=rbind(tmp,tmp_riskscore)

save(tmp,file='Rdata/tcga/uni_cox_clinical.Rdata')

write.csv(tmp,file = 'export/tcga/uni_cox_clinical.csv')
rownames(tmp)
s = as.formula(paste("Surv(time, event) ~ ",paste(c(
  "age",
  "gender",
  "family_history",
  "inflammation_grade",
  "stage",
  "riskscore"),collapse = "+")))
s
model <- coxph(s, data = meta )
summary(model)$concordance
dev.off()
png('plot/tcga/all_clinical_forest.png',height = 400,width = 400)
ggforest(model, data =meta)
dev.off()
save(model,s,file = 'Rdata/tcga/clinical_cox.Rdata')

rm(list=ls())
load('Rdata/tcga/clinical_cox.Rdata')
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

dat2$Trait=c('age',"gender",'family_history',
             "inflammation_grade","stage", 
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
options(digits = 2)
dat2$Coef=c("Coef","0.004","0.214", "0.588","0.262",
            "0.147","0.897")
str(dat2)
dat2$HR=as.numeric(dat2$HR)
dat2$lower=as.numeric(dat2$lower)
dat2$upper=as.numeric(dat2$upper)

colnames(dat2)
png('plot/tcga/multi_clinical_cox.png',height = 300,width =500 )
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
  is.summary = c(T, rep(F, 6)),
  align = "l",
  hrzl_lines = list(
    "1" = gpar(lty=1.5),
    "2" = gpar(lty=1.5),
    "8"= gpar(lty=1.5)),
  colgap = unit(5, 'mm')
)
dev.off()

###单因素
rm(list=ls())
uni_clinical=read.csv('export/tcga/pre/uni_cox_clinical_zhengli.csv')
tmp=uni_clinical

tmp=rbind(c("Trait","Coef" ,NA, NA, NA, "HR", "P"),tmp)
str(tmp)
tmp$HR=as.numeric(tmp$HR)
tmp$lower=as.numeric(tmp$lower)
tmp$upper=as.numeric(tmp$upper)
dev.off()
library(forestplot)
png('plot/tcga/pre/unni_clinical_cox.png',height = 300,width =500 )
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
  is.summary = c(T, rep(F, 6)),
  align = "l",
  hrzl_lines = list(
    "1" = gpar(lty=1.5),
    "2" = gpar(lty=1.5),
    "8"= gpar(lty=1.5)),
  colgap = unit(5, 'mm')
)
dev.off()


rm(list = ls())
load('Rdata/tcga/pre/clinicaldata_juzhenshuju.Rdata')
identical(rownames(meta),colnames(exp_dat))
exp_dat=exp_dat[,match(rownames(meta),colnames(exp_dat))]
identical(rownames(meta),colnames(exp_dat))
Group=meta$ri
Group=factor(Group,levels = c('lowrisk','highrisk'))
pca.plot = draw_pca(exp_dat,Group,addEllipses = F);pca.plot
ggsave('plot/tcga/pre/pca_risk.png',height = 5,width = 5)
