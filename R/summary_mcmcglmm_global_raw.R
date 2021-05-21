library("MCMCglmm")
library("ape")
library("tidyverse")
library("reshape2")
library("cowplot")
library("RColorBrewer")

source("../R/data_preparation_mcmcglmm.R")

load("../output/mcmcglmm_global_pca_raw.RData")

names(fullmcmc.results.global.pca) <- paste0("tree", 1:112)

epsnet.raw.table <-
    bind_rows(lapply(fullmcmc.results.global.pca,
                     function(x){res <- as.data.frame(x$Sol); names(res) = c("Intercept", "Extinction Fraction", "NetDiversification", "PC1", "PC2", "Extinction Fraction:NetDiversification", "Extinction Fraction:PC1", "Extinction Fraction:PC2", "NetDiversification:PC1", "NetDiversification:PC2")
                         return(res)}))

epsnet.raw.table$tree <- rep(gsub("epsilon.", "", names(fulldata.net.bin)[26:137]), each = 1000)

epsnet.raw.plot <- melt(epsnet.raw.table)

posterior.quants.epsnet.raw <- data.frame(
    variable = unique(epsnet.raw.plot$variable),
    aggregate(epsnet.raw.plot$value, by = list(epsnet.raw.plot$variable), FUN = function(x){quantile(x, c(0.025, 0.1, 0.5, 0.9, 0.975))})$x
)
names(posterior.quants.epsnet.raw) <- c("variable", "quant2.5", "quant10", "median", "quant90", "quant97.5")

medians.per.tree.epsnet.raw <- aggregate(epsnet.raw.plot$value, by = list(epsnet.raw.plot$tree, epsnet.raw.plot$variable), FUN = median)
names(medians.per.tree.epsnet.raw) <- c("tree", "variable", "value")
