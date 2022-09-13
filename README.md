# TIGENS

# Installation of TIGENS
#install.packages("devtools")

devtools::install_github("wenchaoli1007/TIGENS")




# Example

library(TIGENS)

library(dplyr)

library(readr)


setwd("~/TIGENS")


sepsis_valerie_degs = read_tsv("../data/covid2_severe_mild_degs.txt")

sepsis_valerie_degs = sepsis_valerie_degs %>% filter(abs(avg_log2FC) > 0.25)


inputs = sepsis_valerie_degs[,c("gene", "avg_log2FC")]

colnames(inputs) = c("gene", "value")

load("./test/TIGs.RData")


refs = data.frame(ti_gene)

res = DisAssoc(inputs, "down", c(1), 0, "no", refs)


plots = PlotAssoc(res)

plots
