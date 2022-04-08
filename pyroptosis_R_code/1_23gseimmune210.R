rm(list=ls())
library(IMvigor210CoreBiologies)
#load('Rdata/lass_cox_model.Rdata')
load('Rdata/tcga/pre/model.Rdata')
data(cds)
#save(cds,file = 'Rdata/tcga/IMvigor210CoreBiologies.Rdata')
counts = counts(cds)
counts[1:4,1:4]
an = fData(cds)
an[1:4,1:4]
pd = pData(cds)
pd[1:4,1:4]
counts = counts[,match(rownames(pd),colnames(counts))]
an = an[!duplicated(rownames(an)),]
counts = counts[!duplicated(rownames(counts)),]

g = intersect(an$entrez_id,rownames(counts))
an = an[g,]
counts = counts[g,]

k = (!duplicated(an$symbol))&(!is.na(an$symbol));table(k)
counts = counts[k,]
an = an[k,]
rownames(counts) = an$symbol
meta = pd[,c("os","censOS","Best Confirmed Overall Response",'Sex')]
colnames(meta)[3] = "Response"
meta$Response[meta$Response=="NE"]=NA

str(meta)
colnames(meta)[1:2] = c("time","event")
  exp = log2(edgeR::cpm(counts)+1)
dim(exp)
exp = exp[gs,]
identical(rownames(meta),colnames(exp))
head(meta)
library(survminer)
gs

rownames(exp)
class(exp)
exp=as.data.frame(t(exp))

riskscore=predict(model,exp)

names(riskscore)=rownames(exp)
identical(rownames(exp),rownames(meta))
meta$riskscore=riskscore

res.cut = surv_cutpoint(meta, time = "time", event = "event", variables = "riskscore",minprop = 0.3)
cut = res.cut[["cutpoint"]][1, 1]
cut
plot(res.cut, "riskscore", palette = "npg")

meta$ri=ifelse(meta$riskscore<cut,'lowrisk','highrisk')
meta$ri=factor(meta$ri,levels = c('lowrisk','highrisk'))
ri=meta$ri
table(ri)
save(meta,file='Rdata/immune210_clinical.Rdata')
library(survival)
library(survminer)
sfit <- survfit(Surv(time, event)~ri, data=meta)
summary(sfit)$table
p_sur_gseim210=ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
                          risk.table =TRUE,pval =TRUE,
                          conf.int =TRUE,xlab ="Time in months", 
                          ggtheme =theme_light(), 
                          ncensor.plot = TRUE)
library(patchwork)
p1=p_sur_gseim210$plot
p2=p_sur_gseim210$table
p3=p_sur_gseim210$ncensor.plot
p1 /p2  + plot_layout(design = 'A
                         A
                         A
                         A
                         B')

ggsave(file='plot/gseim210_kmplot.png',height = 5,width = 5)


library(ggpubr)
table(meta$Response)
dat1 = na.omit(meta)
table(dat1$Response)
ggboxplot(data= dat1,x = "Response",y = "riskscore",color = "Response",add = "jitter")+
  stat_compare_means(label.y = 4)+
  stat_compare_means(comparisons = list(c("CR","PR"),
                                        c("CR","SD"),
                                        c("CR","PD"),
                                        c("PR","SD"),
                                        c("PR","PD"),
                                        c("SD","PD")))
ggsave(file='plot/gseim210_drawplot.png',height = 5,width = 5)
a=dat1[dat1$Response=='CR',]
b=dat1[dat1$Response=='PR',]
e=dat1[dat1$Response=='SD',]
f=dat1[dat1$Response=='PD',]
summary(a$riskscore)
summary(b$riskscore)
summary(e$riskscore)
summary(f$riskscore)

library(dplyr)
meta$Response2 = ifelse(meta$Response %in% c("SD","PD"),"SD/PD","CR/PR")
class(meta$Response2)
meta$Response2=factor(meta$Response2,levels = c("CR/PR","SD/PD"))
meta$ri=factor(meta$ri,levels = c("lowrisk","highrisk"))
class(meta$ri)
table(meta$ri)

p_reponse=ggstatsplot::ggbarstats(
  data = meta,
  x = Response,
  y = ri,
  ggtheme = ggplot2::theme_classic(),
  palette = "category10_d3", 
  package = "ggsci", 
  k = 3, 
  perc.k = 1 ,
  results.subtitle = T,xlab ='',
  return='plot',messages=T,
  bf.message = F
) 
ggsave(file='plot/gseim210_percent.png',height = 5,width = 5)
colnames(meta)
dat_time = meta[,c(1,2,4)]
                 
head(dat_time)
library(survminer)
library(survival)
library(timeROC)
result <-with(dat_time, timeROC(T=time,
                                delta=event,
                                marker=riskscore,
                                cause=1,
                                times=c(12,24),
                                iid = TRUE))
plot_dat = data.frame(fpr = as.numeric(result$FP),
                      tpr = as.numeric(result$TP),
                      time = rep(as.factor(c(12,24)),each = nrow(result$TP)))
library(ggplot2)
ggplot() + 
  geom_line(data = plot_dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582"),
                     labels = paste0("AUC of ",c(1,2),"-y survival: ",
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
ggsave('plot/im210_timeROC.png',height = 5,width = 5)


