library("MCMCglmm")
library("ape")
library("tidyverse")
library("reshape2")
library("cowplot")
library("RColorBrewer")

source("../R/data_preparation_mcmcglmm.R")

load("../output/mcmcglmm_lambdamu_pca_raw.RData")

names(fullmcmc.results.lambdamu.pca) <- paste0("tree", 1:112)

lambdamu.raw.table <-
    bind_rows(lapply(fullmcmc.results.lambdamu.pca,
                     function(x){res <- as.data.frame(x$Sol); names(res) = c("Intercept", "Speciation", "Extinction", "PC1", "PC2", "Speciation:Extinction", "Speciation:PC1", "Speciation:PC2", "Extinction:PC1", "Extinction:PC2")
                         return(res)}))

lambdamu.raw.table$tree <- rep(gsub("epsilon.", "", names(fulldata.net.bin)[26:137]), each = 1000)

lambdamu.raw.plot <- melt(lambdamu.raw.table)
levels(lambdamu.raw.plot$variable) <- levels(lambdamu.raw.plot$variable)[c(1:6, 7:10)]

posterior.quants.lambdamu.raw <- data.frame(
    variable = unique(lambdamu.raw.plot$variable),
    aggregate(lambdamu.raw.plot$value, by = list(lambdamu.raw.plot$variable), FUN = function(x){quantile(x, c(0.025, 0.1, 0.5, 0.9, 0.975))})$x
)
names(posterior.quants.lambdamu.raw) <- c("variable", "quant2.5", "quant10", "median", "quant90", "quant97.5")

medians.per.tree.lambdamu.raw <- aggregate(lambdamu.raw.plot$value, by = list(lambdamu.raw.plot$tree, lambdamu.raw.plot$variable), FUN = median)
names(medians.per.tree.lambdamu.raw) <- c("tree", "variable", "value")
