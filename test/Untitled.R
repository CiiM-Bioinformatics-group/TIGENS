library(TIGENS)
library(dplyr)
library(readr)

setwd("~/TIGENS")

inputs = data.frame(read_tsv("./degs.txt"))

inputs = inputs[,c("gene", "avg_log2FC", "patient")]
load("./test/TIGs.RData")

Disease_Name = "Sepsis"
Sub_Type = "Immune paralysis"
Reg_Dir = "up"
Test_Dis_Cutoff = 0.25
Tested_TI_Cluster = 4
Tested_TI_Cluster2 = "TM4"
Pval_Cutoff = 0.1

refs = data.frame(ti_gene)
res = DisAssoc(disease_data = inputs,
               reg_dir = Reg_Dir,
               cluster_id = Tested_TI_Cluster,
               cutoff = Test_Dis_Cutoff,
               sub_type = Sub_Type,
               ti_gene = refs,
               pval_cutoff = Pval_Cutoff)

## basic version
plots = GSEAplot_Basic(res, 1,
                       disease_name = Disease_Name,
                       sub_type = Sub_Type,
                       cluster_id = Tested_TI_Cluster2,
                       reg_dir = Reg_Dir)
plots


## 
plots2 = GSEAplot(res, 1,
                   disease_name = Disease_Name,
                   sub_type = Sub_Type,
                   cluster_id = Tested_TI_Cluster2,
                   reg_dir = Reg_Dir,
                   pvalue_table = TRUE)

plots2



genes_to_enrich = enriched_gene_to_pathway(res)

data.frame(res@result)


