rm(list = ls())
a=read.csv('export/tcga/pre/unicox_tcga_zhengli.csv')
tmp=a

tmp=rbind(c("Trait","Coef" ,NA, NA, NA, "HR", "P"),tmp)
str(tmp)
tmp$HR=as.numeric(tmp$HR)
tmp$lower=as.numeric(tmp$lower)
tmp$upper=as.numeric(tmp$upper)
dev.off()
library(forestplot)
png('plot/tcga/pre/tcga_unni_cox.png',height = 300,width =500 )
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
  is.summary = c(T, rep(F, 11)),
  align = "l",
  hrzl_lines = list(
    "1" = gpar(lty=1.5),
    "2" = gpar(lty=1.5),
    "13"= gpar(lty=1.5)),
  colgap = unit(5, 'mm')
)
dev.off()


