library(TIGENS)
library(dplyr)
library(readr)

setwd("~/hzi-li/project/300BCG/git/TIGENS")

inputs = data.frame(read_tsv("../../algorithm/results/disease/sepsis_valerie_degs.txt"))

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


plots = GSEAplot_Basic(res, 1,
                       disease_name = Disease_Name,
                       sub_type = Sub_Type,
                       cluster_id = Tested_TI_Cluster2,
                       reg_dir = Reg_Dir)
plots


plots2 = GSEAplot(res, 1,
                   disease_name = Disease_Name,
                   sub_type = Sub_Type,
                   cluster_id = Tested_TI_Cluster2,
                   reg_dir = Reg_Dir,
                   pvalue_table = TRUE)

plots2



genes_to_enrich = enriched_gene_to_pathway(res)

data.frame(res@result)


##### test
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList

  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]

  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  if (length(object@gene2Symbol) == 0) {
    df$gene <- names(geneList)
  } else {
    df$gene <- object@gene2Symbol[names(geneList)]
  }

  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}


a = gsInfo(res, geneSetID = 1)



gseaplot(res, geneSetID = 1) +
  theme(axis.title.y.left = element_text(size = 0.8)) +
  geom_gsea_gene(genes_to_enrich,a)


