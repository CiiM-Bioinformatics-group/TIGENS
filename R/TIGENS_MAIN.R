#' Pre-Processing of input dataset
#'
#' input could be interested gene list. Each row is one gene and columns are "gene" and "value".
#' "value" could be log2FC or specific value.
#' colnames of other metadata should be "subtype".
#' @param disease_data data.frame with colnames "gene", "value", ("subtype")
#' @param sub_type which disease subtype should be tested
#' @import dplyr
#' @importFrom dplyr %>%
#' @return Return ranked data.frame based on "value"
#'
#' @export

PreProDisData <- function(disease_data, sub_type)
{
  if (ncol(disease_data) == 2)
    disease_data$subtype = "all"
  else if(ncol(disease_data) > 2)
    disease_data = disease_data %>% filter(subtype == sub_type)
  else if(ncol(disease_data) < 2)
    stop("error: at least two columns (gene, value) should be provided!")

  rank_dd = disease_data[order(disease_data$value, decreasing = TRUE), ]

  return(rank_dd)

}


#' Check the association between specific disease and monocytes sub-populations
#'
#' This tool will use GASE embedded in clusterProfiler to check the association between a specific
#' disease and identified monocytes sub-populations. For each input dataset, diseases could be
#' classified into sub-types. For each disease/sub-type, genes will be ranked based on lof2FC.
#' 4 sub-populations would be calculated one by one.
#' @param disease_data The input data.frame with colnames "gene", "value", ("subtype")
#' @param reg_dir Value should be "up" or "down"
#' @param cluster_id  Numeric variable or vector
#' @param cutoff  Select genes to gsea
#' @param sub_type which disease subtype should be tested
#' @param ti_gene reference database downloaded from github
#' @import clusterProfiler
#' @import dplyr
#' @importFrom clusterProfiler GSEA
#' @importFrom dplyr %>%
#' @return Enriched GSEA results for plotting or other downstream analysis

#'
#' @export


DisAssoc = function(disease_data, reg_dir, cluster_id, cutoff, sub_type, ti_gene)
{

  rank_dd = PreProDisData(disease_data, sub_type)
  if (reg_dir == "up")
    dir_ti_gene = ti_gene %>% filter(avg_log2FC > 0 & cluster == cluster_id)

  else if(reg_dir == "down")
    dir_ti_gene = ti_gene %>% filter(avg_log2FC < 0 & cluster == cluster_id)
  else
    stop("error: please add regulatory dirction!")


  gsea_results = calgsea(dir_ti_gene, cluster_id, rank_dd, rank_dd$subtype)


  return(gsea_results)

}



#' Calcualte GSEA enrichment
#' @param dir_ti_gene Data.frame includes TIGs, avg_lo2FC, regulated direction
#' @param j The tested sub-population
#' @param rank_dd Ranked disease DEGs
#' @param sub_type The tested sub-type of disease
#' @import clusterProfiler
#' @importFrom clusterProfiler GSEA
#' @return GSEA enrichment results
#'
#' @export


calgsea = function(dir_ti_gene, j, rank_dd, sub_type)
{
  term_gene = data.frame(ti_cluster = j,
                         geneID = dir_ti_gene$gene)
  rank_dd = rank_dd %>% filter(subtype == sub_type)
  ranks_gene = rank_dd$value
  names(ranks_gene) = rank_dd$gene
  gsea_results = GSEA(geneList = ranks_gene,
                      pvalueCutoff = 0.1,
                      minGSSize = 3,
                      TERM2GENE = term_gene,
                      eps = 0)

  return(gsea_results)

}


#' Plot enrichment results
#'
#' @param gsea_results The result from GSEA enrichment
#' @import clusterProfiler
#' @import ggplot2
#' @importFrom clusterProfiler gseaplot
#' @importFrom ggplot2 theme
#' @return gseaplots
#' @export


PlotAssoc = function(gsea_results)
{
  gseaplots = gseaplot(gsea_results, geneSetID = 1) +
    theme(axis.title.y.left = element_text(size = 1))
  return(gseaplots)
}



#' export enriched genes for downstream analysis
#' @param gsea_results The results from GSEA enrichment
#' @export

enriched_gene_to_pathway = function(gsea_results)
{
  genes_to_enrich = data.frame(gsea_results@result)$core_enrichment
  genes_to_enrich = data.frame(genes = strsplit(genes_to_enrich, split = "/")[[1]])

  return(genes_to_enrich)
}

