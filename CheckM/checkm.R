library(ggplot2)
library(foreach)
library(plyr)
library(reshape2)
source("~/software/metabat/benchmark.R")
res <- list(maxbin2=calcPerfBySCG("../checkm_maxbin.txt"), 
			graphbin=calcPerfBySCG("edge_graph_graphbin.txt"), 
			amgl_remove=calcPerfBySCG("final.txt"), 
			maxbin2_magpurify=calcPerfBySCG("../checkm_maxbin_cleaned.txt"), 
			graphbin_magpurify=calcPerfBySCG("edge_graph_graphbin_magpurify.txt"), 
			amgl_magpurify_remove = calcPerfBySCG("final_magpurify.txt"),
			amgl_magpurify_cleaned=calcPerfBySCG("final_noremove_magpurify_cleaned.txt"))
printPerf(res)
pdf("Performance_By_SCG_All.pdf", width=12, height=4)
plotPerf3(res)
dev.off()