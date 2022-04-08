rm(list = ls())
library(GSEABase)
gene_msig=getGmt('import/geneset_jiaowang.gmt')
gene_msig=gene_msig[["REACTOME_PYROPTOSIS"]]@geneIds
gene_cards=read.csv('import/Pyroptosis_GeneCards-SearchResults.csv')
gene_cards=gene_cards$Gene.Symbol[gene_cards$Relevance.score>7]
gene_file=readxl::read_xlsx('import/pyroptosis.xlsx')
gene_file=gene_file$Genes
pyroptosis=union(gene_msig,union(gene_cards,gene_file))
#pyroptosis=gene_file
table(duplicated(pyroptosis))
table_pyroptosis=data.frame(gene.symbols=pyroptosis)
write.csv(pyroptosis,'export/table_pyroptosis.csv')
save(pyroptosis,file = 'Rdata/pyroptosis.Rdata')
