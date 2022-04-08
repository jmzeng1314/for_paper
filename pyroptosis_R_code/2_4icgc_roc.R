rm(list = ls())
load('Rdata/icgc/pre/linchuangzuiquan_icgc.Rdata')
library(timeROC)
library(survival)
str(meta)
meta$gender=ifelse(meta$gender=='male',1,2)
df2=meta
colnames(df2)
str(df2)
# riskScore的ROC曲线
ROC.risk <- timeROC(T=df2$time,
                    delta=df2$event,   
                    marker=df2$riskscore,   
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


###stage
Roc.stage=timeROC(T=df2$time,   
                  delta=df2$event,   
                  marker=df2$stage,   
                  cause=1,   
                  weighting="marginal",   
                  times=36,   
                  iid=TRUE)

png('plot/icgc/pre/muti_index_roc_all.png',height = 500,width = 500)
plot(ROC.risk, time = 36, col="#E41A1C", lwd=2, title = "")
plot(ROC.gender, time = 36, col="#A65628", lwd=2, add = T)
plot(ROC.age, time = 36, col="#4DAF4A", lwd=2, add = T)
#plot(ROC.inflammation, time = 36, col="#377EB8", lwd=2, add = T)
plot(Roc.stage, time = 36, col="#984EA3", lwd=2, add = T)
legend("bottomright",
       c(paste0("riskscore:",round(ROC.risk[["AUC"]][2],2)), 
         paste0("gender:",round(ROC.gender[["AUC"]][2],2)), 
         paste0("age:",round(ROC.age[["AUC"]][2],2)),
        # paste0("inflammation_grade:",round(ROC.inflammation[["AUC"]][2],2)),
         paste0("stage:",round(Roc.stage[["AUC"]][2],2))
         
       ),
       col=c("#E41A1C", "#A65628", "#4DAF4A","#984EA3"),
       lty=1, lwd=2,bty = "n")
dev.off()

ROC <- timeROC(T=df2$time,   
               delta=df2$event,   
               marker=df2$riskscore,   
               cause=1,                #阳性结局指标数值
               weighting="marginal",   #计算方法，默认为marginal
               times=c(12, 24, 36),       #时间点，选取1年，3年和5年的生存率
               iid=TRUE)
png('plot/icgc/pre/timeROC_all.png',height = 500,width = 500)
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


