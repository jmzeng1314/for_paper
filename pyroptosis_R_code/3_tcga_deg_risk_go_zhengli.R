rm(list = ls())
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
library(ggthemes)
load("Rdata/tcga/pre/deg_risk_GO.Rdata")
a=c('GO:0006874','GO:0071248','GO:0055074','GO:0045165','GO:0007204',
    'GO:1990351','GO:1902495','GO:0062023','GO:0016324','GO:0031225',
    'GO:0005201','GO:0008528','GO:0016493','GO:0015108','GO:0019957')
ego@result=ego@result[ego@result$ID%in%a,]
dim(ego@result)
barplot(ego)
dev.off()

png('plot/tcga/pre/deg_risk_tcga_GO.png',height = 800,width = 800)
dotplot(ego, split = "ONTOLOGY", font.size = 10, 
        showCategory = 5) + 
  facet_grid(ONTOLOGY ~ ., scale = "free") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45))
dev.off()
