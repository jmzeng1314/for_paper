rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
proj="TCGA-LIHC"
if(!dir.exists("clinical"))dir.create("clinical")
if(!dir.exists("expdata"))dir.create("expdata")
dir()
length(dir("./clinical/"))
length(dir("./expdata/"))

library(XML)
result <- xmlParse("./clinical/00a9a7f4-06eb-40fb-8d3a-66f5f5d315f7/nationwidechildrens.org_clinical.TCGA-KR-A7K0.xml")
rootnode <- xmlRoot(result)
xmlSize(rootnode)
xmldataframe <- xmlToDataFrame(rootnode[2])
head(t(xmlToDataFrame(rootnode[2])))
xmls = dir("clinical/",pattern = "*.xml$",recursive = T)
cl = list()
for(i in 1:length(xmls)){
  result <- xmlParse(paste0("clinical/",xmls[[i]]))
  rootnode <- xmlRoot(result)
  cl[[i]] = xmlToDataFrame(rootnode[2])
}
clinical <- do.call(rbind,cl)
clinical[1:3,1:3]


options(stringsAsFactors = F)
x = read.table("expdata/00eb7cff-2a7d-464e-92f2-04b796a79181/3243c347-c219-487f-87fa-ec4046feee81.htseq.counts.gz",
               row.names = 1,sep = "\t")
x2 = read.table("expdata/0a4f6166-47ea-4c16-9fb1-e8276cb2a9cb/7ab845f8-320e-4c93-9d5c-bf66bb8d121c.htseq.counts.gz",
                row.names = 1,sep = "\t")
identical(rownames(x),rownames(x2))

count_files = dir("expdata/",pattern = "*.htseq.counts.gz$",recursive = T)
exp = list()
for(i in 1:length(count_files)){
  exp[[i]] <- read.table(paste0("expdata/",count_files[[i]]),row.names = 1,sep = "\t")
}
exp <- do.call(cbind,exp)
dim(exp)
exp[1:4,1:4]



meta <- jsonlite::fromJSON("metadata.cart.2021-08-21.json")
colnames(meta)
class(meta$associated_entities)
meta$associated_entities[[1]]
meta$associated_entities[[1]]$entity_submitter_id

ID = sapply(meta$associated_entities,
            function(x){x$entity_submitter_id})
file2id = data.frame(file_name = meta$file_name,
                     ID = ID)

head(file2id$file_name,2)
head(count_files,2)
count_files2 = stringr::str_split(count_files,"/",simplify = T)[,2]
table(count_files2 %in% file2id$file_name)

file2id = file2id[match(count_files2,file2id$file_name),]
colnames(exp) = file2id$ID
#表达矩阵整理完成,需要过滤一下那些在很多样本里表达量都为0的基因。过滤标准不唯一。
dim(exp)
exp = exp[apply(exp, 1, function(x) sum(x > 1) > 0.25*ncol(exp)), ]#筛选在10个样本里表达量大于1的样本
dim(exp)
exp[1:4,1:4]
exp = as.matrix(exp)
class(exp[,1])


### 分组信息
#根据样本ID的第14-15位，给样本分组（tumor和normal）
table(str_sub(colnames(exp),14,15))
#02是血液样本的意思，所以要把血液样本去掉
table(str_sub(colnames(exp),14,15)!='02')
exp=exp[,str_sub(colnames(exp),14,15)!='02']
Group = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')
Group = factor(Group,levels = c("normal","tumor"))
table(Group)
save(exp,clinical,Group,proj,file = paste0(proj,"_gdc.Rdata"))














