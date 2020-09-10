#!/usr/bin/env Rscript

# install.packages("devtools")
# devtools::install_github("hemberg-lab/scfind")

library(scfind)
library(SingleCellExperiment)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(UpSetR)
library(ComplexHeatmap)
library(pheatmap)
library(cowplot)


# standardize theme
main <- theme_minimal() + theme(axis.title.y = element_text(size=12), axis.title.x = element_text(size=12), text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) 
colorcode = "Dark2"
RColorBrewer::brewer.pal(8,colorcode)
dotsize = .7
linesize = .3
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
tritPalette <- c("#519FAC", "#597C86", "#C84B48","#d8b365", "#f5f5f5", "#5ab4ac", "#EE6867", "#F4AAB2", "#B48B25", "#807271")
colour18 <- unique(c(brewer.pal(9, "Set3"), brewer.pal(9, "Paired")))
colourGO <- c(cbPalette, tritPalette, "black")

legend_order <- c("BCA", "MCA", "TM, 10X", "TM, FACS", "sciATACseq", "MOCA")

scale_y_log <- scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
