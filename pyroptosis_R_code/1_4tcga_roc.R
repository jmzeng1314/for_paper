rm(list=ls())
load('Rdata/tcga/pre/clinicaldata_juzhenshuju.Rdata')
library(timeROC)
library(survival)
str(meta)
meta$gender=ifelse(meta$gender=='MALE',1,2)
meta$stage=ifelse(meta$stage=='I',1,
                  ifelse(meta$stage=='II',2,
                         ifelse(meta$stage=='III',3,
                                ifelse(meta$stage=='IV',4,NA))))

meta$family_history=ifelse(meta$family_history=='NO',0,
                           ifelse(meta$family_history=='YES',1,NA))

meta$inflammation_grade=ifelse(meta$inflammation_grade=='G1',1,
                               ifelse(meta$inflammation_grade=='G2',2,
                                      ifelse(meta$inflammation_grade=='G3',3,
                                             ifelse(meta$inflammation_grade=='G4',4,NA))))                          
meta$adjacent_inflammation=ifelse(meta$adjacent_inflammation=='None',0,
                                  ifelse(meta$adjacent_inflammation=='Mild',1,
                                         ifelse(meta$adjacent_inflammation=='Severe',2,NA)))
colnames(meta)
meta=meta[,c(1,2,6,7,8,9,11,17,20,22,23,24,25)]
str(meta)

df2=meta
colnames(df2)=c('ID','event','age','gender','stage','family_history',
                'inflammation_grade','adjacent_inflammation',
                'time','T','N','M','risk_score')
str(df2)
# riskScore的ROC曲线
ROC.risk <- timeROC(T=df2$time,
                    delta=df2$event,   
                    marker=df2$risk_score,   
                    cause=1,                
                    weighting="marginal",   
                    times=36,   
                    iid=TRUE)

# gender的ROC曲线
ROC.gender <- timeROC(T=df2$time,   
                      delta=df2$event,   
                      marker=df2$gender,   
                      cause=1,   
                      weighting="marginal",   
                      times=36,   
                      iid=TRUE)
# age的ROC曲线
ROC.age <- timeROC(T=df2$time,   
                   delta=df2$event,   
                   marker=df2$age,   
                   cause=1,   
                   weighting="marginal",   
                   times=36,   
                   iid=TRUE)

##inflammation_grade
ROC.inflammation=timeROC(T=df2$time,   
                         delta=df2$event,   
                         marker=df2$inflammation_grade,   
                         cause=1,   
                         weighting="marginal",   
                         times=36,   
                         iid=TRUE)
###stage
Roc.stage=timeROC(T=df2$time,   
                  delta=df2$event,   
                  marker=df2$stage,   
                  cause=1,   
                  weighting="marginal",   
                  times=36,   
                  iid=TRUE)

png('plot/tcga/pre/muti_index_roc_all.png',height = 500,width = 500)
plot(ROC.risk, time = 36, col="#E41A1C", lwd=2, title = "")
plot(ROC.gender, time = 36, col="#A65628", lwd=2, add = T)
plot(ROC.age, time = 36, col="#4DAF4A", lwd=2, add = T)
plot(ROC.inflammation, time = 36, col="#377EB8", lwd=2, add = T)
plot(Roc.stage, time = 36, col="#984EA3", lwd=2, add = T)
legend("bottomright",
       c(paste0("riskscore:",round(ROC.risk[["AUC"]][2],2)), 
         paste0("gender:",round(ROC.gender[["AUC"]][2],2)), 
         paste0("age:",round(ROC.age[["AUC"]][2],2)),
         paste0("inflammation_grade:",round(ROC.inflammation[["AUC"]][2],2)),
         paste0("stage:",round(Roc.stage[["AUC"]][2],2))
         
       ),
       col=c("#E41A1C", "#A65628", "#4DAF4A","#377EB8","#984EA3"),
       lty=1, lwd=2,bty = "n")
dev.off()

ROC <- timeROC(T=df2$time,   
               delta=df2$event,   
               marker=df2$risk_score,   
               cause=1,                #阳性结局指标数值
               weighting="marginal",   #计算方法，默认为marginal
               times=c(12, 24, 36),       #时间点，选取1年，3年和5年的生存率
               iid=TRUE)
png('plot/tcga/pre/timeROC_all.png',height = 500,width = 500)
plot(ROC, 
     time=12, col="#92C5DE", lwd=2, title = "")   #time是时间点，col是线条颜色
plot(ROC,
     time=24, col="#F4A582", add=TRUE, lwd=2)    #add指是否添加在上一张图中
plot(ROC,
     time=36, col="#66C2A5", add=TRUE, lwd=2)
#添加标签信息
legend("bottomright",
       c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
         paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
         paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2))),
       col=c("#92C5DE", "#F4A582", "#66C2A5"),
       lty=1, lwd=2,bty = "n")
dev.off()


