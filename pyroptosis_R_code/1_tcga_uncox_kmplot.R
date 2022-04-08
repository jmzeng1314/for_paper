rm(list = ls())
a=read.csv('export/tcga/pre/unicox_tcga_zhengli.csv')
colnames(a) = c('Trait','Coef',"HR", "lower", "upper",'HR2','p')
dat2 = rbind(
  c("Trait","Coef" ,NA, NA, NA, "HR", "P"),
  a[1, ],
  a[2:nrow(a), ]
)
str(dat2)
dat2$HR=as.numeric(dat2$HR)
dat2$lower=as.numeric(dat2$lower)
dat2$upper=as.numeric(dat2$upper)

colnames(dat2)
tmp=dat2
#BiocManager::install('forestplot')
library(forestplot)
png('plot/tcga/pre/uni_cox.png',height = 300,width =500 )
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
  is.summary = c(T, rep(F, 3)),
  align = "l",
  hrzl_lines = list(
    "1" = gpar(lty=1.5),
    "2" = gpar(lty=1.5),
    "5"= gpar(lty=1.5)),
  colgap = unit(5, 'mm')
)
dev.off()

