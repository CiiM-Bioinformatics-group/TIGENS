# TIGENS

Test the enrichment of gene signatures from different trained immunity monocyte sub-populations in patients’ transcriptome data.

# Installation of TIGENS
#install.packages("devtools")

devtools::install_github("CiiM-Bioinformatics-group/TIGENS")




# Introduction to TIGENS:

We investigated the heterogeneity of TI induction by single-cell RNA-seq of immune cells collected from 39 healthy individuals before and three months after BCG vaccination, and after ex vivo heterologous stimulation with lipopolysaccharide (LPS). We observed that not only monocytes but also CD8+ T cells showed heterologous transcriptional responses after stimulation with LPS, with an active crosstalk between these two cell types. Unsupervised clustering identified four distinct functional monocyte sub-populations. Enrichment analysis indicated that the interferon-gamma pathway was crucial in BCG-induced TI, and this pathway was up-regulated in functional high-responders. To investigate the characteristics of newly defined TI sub-populations, data-driven and public dataset-based analyses were combined with functional experiments, and revealed STAT1 to be an enriched transcription factor for TI shared in all sub-populations. Finally, we reported the role of type I interferon-related and neutrophil-related TI transcriptional programs in patients with sepsis. These findings provided comprehensive insights into the importance of monocyte heterogeneity during TI in humans. The association of these TI programs with disease could be informative for novel immunotherapy design. 

We use the single-cell RNAseq data generated from Sepsis patients[1], with immune paralysis and macrophage activation-like syndrome (MALS) profiles.

# Below is an example of the enrichment:

We use the DEGs between immune paralysis patients and healthy people and test the up-regulated TIGs in TM4.



# How to use:

The example could be found here: /test/Untitled.R

[1] Valerie Koeken, Inge Grondman, Athanasios Karageorgos, Wenchao Li, Nikolaos Antonakos, Bowen Zhang, Georgia Damoraki, Cheng-Jian Xu, Evangelos J. Giamarellos-Bourboulis, Yang Li, Mihai G. Netea, Single-cell transcriptomics differentiates hyperinflammation from immune paralysis in sepsis patients. (submitted) (2022).
