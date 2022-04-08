rm(list = ls())
load('Rdata/tcga/pre/model.Rdata')
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
dat2$Coef=c("Coef","0.357","0.397", "-0.364","0.301")
str(dat2)
dat2$HR=as.numeric(dat2$HR)
dat2$lower=as.numeric(dat2$lower)
dat2$upper=as.numeric(dat2$upper)

colnames(dat2)
png('plot/tcga/multi_cox_4genes.png',height = 300,width =500 )
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
