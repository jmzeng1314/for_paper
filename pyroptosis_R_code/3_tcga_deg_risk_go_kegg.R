rm(list = ls())
load('Rdata/tcga/pre/tcga_deg_risk.Rdata')
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
library(ggthemes)
deg=rownames(DEG1)[DEG1$change!='NOT']
gene.df <- bitr(deg, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                toType = c( "ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Hs.eg.db)
head(gene.df)
save(gene.df,file='Rdata/tcga/pre/deg_risk_gene_entrezid.Rdata')
table(deg%in%gene.df$SYMBOL)
dim(gene.df)
deg_deg1=DEG1[rownames(DEG1)%in%gene.df$SYMBOL,]
dim(gene.df)
gene.df=gene.df[!duplicated(gene.df$SYMBOL),]
dim(deg_deg1)
identical(rownames(deg_deg1),gene.df$SYMBOL)
deg_deg1$ENTREZID= gene.df$ENTREZID
# 1.GO 富集分析----

#(1)输入数据
deg=deg_deg1
gene_up = deg[deg$change == 'UP','ENTREZID'] 
gene_down = deg[deg$change == 'DOWN','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
#(2)富集
#以下步骤耗时很长，设置了存在即跳过
if(!file.exists("Rdata/tcga/pre/deg_risk_GO.Rdata")){
  ego <- enrichGO(gene = gene_diff,
                  OrgDb= org.Hs.eg.db,
                  ont = "ALL",
                  readable = TRUE)
  ego_BP <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "BP",
                     readable = TRUE)
  ego_MF <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "MF",
                     readable = TRUE)
  ego_CC <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "cc",
                     readable = TRUE)
  #ont参数：One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
  save(ego,ego_BP,ego_MF,ego_CC,file = "Rdata/tcga/pre/deg_risk_GO.Rdata")
}
load("Rdata/tcga/pre/deg_risk_GO.Rdata")





png('plot/tcga/pre/deg_risk_tcga_GO.png',height = 800,width = 800)
dotplot(ego, split = "ONTOLOGY", font.size = 10, 
        showCategory = 5) + 
  facet_grid(ONTOLOGY ~ ., scale = "free") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45))
dev.off()
#geneList 用于设置下面图的颜色
geneList = deg$log2FoldChange
names(geneList)=deg$ENTREZID

#(3)展示top通路的共同基因，要放大看。
#Gene-Concept Network
# cnetplot(ego,categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
# png('plot/deg/Go_deg_patheay.png',height = 600,width =1000)
# cnetplot(ego,foldChange=geneList, circular = TRUE, colorEdge = TRUE)
# dev.off()
#Enrichment Map,这个函数最近更新过，版本不同代码会不同


# 2.KEGG pathway analysis----
#上调、下调、差异、所有基因
#（1）输入数据
gene_up = deg[deg$change == 'UP','ENTREZID'] 
gene_down = deg[deg$change == 'DOWN','ENTREZID'] 
gene_diff = c(gene_up,gene_down)

#（2）对上调/下调/所有差异基因进行富集分析

if(!file.exists("Rdata/tcga/pre/deg_KEGG.Rdata")){
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa')
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa')
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa')
  save(kk.diff,kk.down,kk.up,file = "Rdata/tcga/pre/deg_KEGG.Rdata")
}
load("Rdata/tcga/pre/deg_KEGG.Rdata")
#(3)看看富集到了吗？https://mp.weixin.qq.com/s/NglawJgVgrMJ0QfD-YRBQg
table(kk.diff@result$pvalue<0.05)
table(kk.up@result$pvalue<0.05)
table(kk.down@result$pvalue<0.05)

#(4)双向图
# 富集分析所有图表默认都是用p.adjust,富集不到可以退而求其次用p值，在文中说明即可

down_kegg <- kk.down@result %>%
  filter(pvalue<0.05) %>% #筛选行
  mutate(group=-1) #新增列

up_kegg <- kk.up@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=1)


library(ggthemes)
library(ggplot2)
up.data=up_kegg
down.data=down_kegg 
dat=rbind(up.data,down.data)
colnames(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 

dat=dat[order(dat$pvalue,decreasing = F),]
dat=dat[!duplicated(dat$Description),]
write.csv(dat,file = 'export/tcga/pre/tcga_kegg.csv')
save(dat,file = 'Rdata/tcga/pre/tcga_kegg.Rdata')

dat_up=dat[2:11,]
dat_down=dat[(nrow(dat)):(nrow(dat)-9),]
dat_small=rbind(dat_up,dat_down)

g_kegg=gk_plot <- ggplot(dat_small,aes(reorder(Description, pvalue), y=pvalue)) +
  geom_bar(aes(fill=factor((pvalue>0)+1)),stat="identity", width=0.7, position=position_dodge(0.7)) +
  coord_flip() +
  scale_fill_manual(values=c("#0072B2", "#B20072"), guide="none") +
  labs(x="", y="" ) +
  theme_pander()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line.x = element_line(size = 0.3, colour = "black"),#x轴连线
        axis.ticks.length.x = unit(-0.20, "cm"),#修改x轴刻度的高度，负号表示向上
        #axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),##线与数字不要重叠
        axis.ticks.x = element_line(colour = "black",size = 0.3) ,#修改x轴刻度的线                         
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust=0),
        panel.background = element_rect(fill=NULL, colour = 'white')
  )


g_kegg
#g_kegg +scale_y_continuous(labels = c(10,5,0,5,10))#这里c()里面的数据要根据数据修改，注意一定要看清楚
ggsave(g_kegg,filename = 'plot/tcga/pre/tcga_kegg_up_down.png',height = 8,width = 8)

write.csv(ego@result,file = 'export/tcga/pre/tcga_deg_go.csv')




