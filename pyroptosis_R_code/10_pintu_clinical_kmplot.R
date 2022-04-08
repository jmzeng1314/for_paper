###boxplot_tcga_icgc
rm(list = ls())
load('Rdata/tcga/pre/tcga_boxplot_clinical.Rdata')
age_tcga=b_age
gender_tcga=b_gender
famhis_tcga=b_familyhistory
stage_tcga=b_stage
inf_tcga=b_inflammation
load('Rdata/icgc/pre/icgc_boxplot_clinical.Rdata')
age_icgc=b_age
gender_icgc=b_gender
stage_icgc=b_stage
library(patchwork)
(age_tcga|age_icgc)/
  (gender_tcga|gender_icgc)/
  (famhis_tcga|inf_tcga)/
  (stage_tcga|stage_icgc)
ggsave('plot/tcga_icgc_clinical_boxplot.png',height = 16,width = 10)


load('Rdata/tcga/pre/train_tcga_km_clinical.Rdata')
p_younger_tcga=p_younger
p_older_tcga=p_older
p_mal_tcga=p_mal
p_femal_tcga=p_femal
p_hist_tcga=p_hist
p_nonhist_tcga=p_nonhist
p_mild_inflammation_tcga=p_mild_inflammation
p_sever_inflammation_tcga=p_sever_inflammation
p_early_stage_tcga=p_early_stage
p_later_stage_tcga=p_later_stage

load('Rdata/icgc/pre/icgc_km_clinical.Rdata')
p_younger_icgc=p_younger
p_older_icgc=p_older
p_mal_icgc=p_mal
p_femal_icgc=p_femal
p_early_stage_icgc=p_early_stage
p_later_stage_icgc=p_later_stage
library(patchwork)
(inf_tcga|stage_tcga|stage_icgc)/
  (p_mild_inflammation_tcga$plot|p_sever_inflammation_tcga$plot|p_early_stage_tcga$plot)/
  (p_later_stage_tcga$plot|p_early_stage_icgc$plot|p_later_stage_icgc$plot)
ggsave('plot/tcga_icgc_clinical_km_boxplot.png',height =15,width = 15 )

(p_younger_tcga$plot|p_younger_icgc$plot|p_older_tcga$plot|p_older_icgc$plot)/
  (p_mal_tcga$plot|p_mal_icgc$plot|p_femal_tcga$plot|p_femal_icgc$plot)/
  (p_hist_tcga$plot|p_nonhist_tcga$plot|p_mild_inflammation_tcga$plot|p_sever_inflammation_tcga$plot)/
  (p_early_stage_tcga$plot|p_early_stage_icgc$plot|p_later_stage_tcga$plot|p_later_stage_icgc$plot)
ggsave('plot/tcga_icgc_clinical_km.png',height = 15,width = 15)

dev.off()  
inf_tcga
ggsave('plot/tcga/pre/inf_box_tcga.png',height = 5,width = 5)
dev.off()
stage_tcga
ggsave('plot/tcga/pre/stage_box_tcga.png',height = 5,width = 5)
dev.off()
stage_icgc
ggsave('plot/icgc/pre/stage_box_icgc.png',height = 5,width = 5)
dev.off()
p_mild_inflammation_tcga$plot
ggsave('plot/tcga/pre/mile_inf_km_tcga.png',height = 5,width = 5)
dev.off()
p_sever_inflammation_tcga$plot
ggsave('plot/tcga/pre/sever_inf_km_tcga.png',height = 5,width = 5)
dev.off()
p_early_stage_tcga$plot
ggsave('plot/tcga/pre/early_stage_km_tcga.png',height = 5,width = 5)
dev.off()
p_later_stage_tcga$plot
ggsave('plot/tcga/pre/later_stage_km_tcga.png',height = 5,width = 5)
dev.off()
p_early_stage_icgc$plot
ggsave('plot/icgc/pre/early_stage_km_icgc.png',height = 5,width = 5)
dev.off()
p_later_stage_icgc$plot
ggsave('plot/icgc/pre/later_stage_km_icgc.png',height = 5,width = 5)
dev.off()
