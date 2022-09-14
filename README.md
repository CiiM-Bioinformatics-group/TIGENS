# TIGENS

Test the enrichment of gene signatures from different trained immunity monocyte sub-populations in patients’ transcriptome data.

# Installation of TIGENS
#install.packages("devtools")

devtools::install_github("CiiM-Bioinformatics-group/TIGENS")




# Introduction to TIGENS:

We investigated the heterogeneity of TI induction by single-cell RNA-seq of immune cells collected from 39 healthy individuals before and three months after BCG vaccination, and after ex vivo heterologous stimulation with lipopolysaccharide (LPS). We observed that not only monocytes but also CD8+ T cells showed heterologous transcriptional responses after stimulation with LPS, with an active crosstalk between these two cell types. Unsupervised clustering identified four distinct functional monocyte sub-populations. Enrichment analysis indicated that the interferon-gamma pathway was crucial in BCG-induced TI, and this pathway was up-regulated in functional high-responders. To investigate the characteristics of newly defined TI sub-populations, data-driven and public dataset-based analyses were combined with functional experiments, and revealed STAT1 to be an enriched transcription factor for TI shared in all sub-populations. Finally, we reported the role of type I interferon-related and neutrophil-related TI transcriptional programs in patients with sepsis. These findings provided comprehensive insights into the importance of monocyte heterogeneity during TI in humans. The association of these TI programs with disease could be informative for novel immunotherapy design. 

The visualization is based on the enricheplot R package developed by Guangchuang YU[1].

We use the single-cell RNAseq data generated from Sepsis patients[2], with immune paralysis and macrophage activation-like syndrome (MALS) profiles.

# Below is an example of the enrichment:

We use the DEGs between immune paralysis patients and healthy people and test the up-regulated TIGs in TM4.

![GSEAexample1](https://user-images.githubusercontent.com/107392830/189920977-233fbbe7-1797-48df-830d-9b96bacc00d9.jpg)

![GSEAexample2](https://user-images.githubusercontent.com/107392830/189921010-c8d1480f-e6a3-4c53-b74a-2995c2495351.jpg)


# How to use:

The example could be found here: /test/Untitled.R

# Reference: 
[1] T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141. doi: 10.1016/j.xinn.2021.100141
[2] Valerie Koeken, Inge Grondman, Athanasios Karageorgos, Wenchao Li, Nikolaos Antonakos, Bowen Zhang, Georgia Damoraki, Cheng-Jian Xu, Evangelos J. Giamarellos-Bourboulis, Yang Li, Mihai G. Netea, Single-cell transcriptomics differentiates hyperinflammation from immune paralysis in sepsis patients. (submitted) (2022).

# Citation:

A single cell view on host immune transcriptional response to in vivo BCG-induced trained immunity. (submitted). 2022.

