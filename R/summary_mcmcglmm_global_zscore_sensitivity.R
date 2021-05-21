library("MCMCglmm")
library("ape")
library("tidyverse")
library("reshape2")
library("cowplot")
library("RColorBrewer")

### Global + PCA

load("../output/mcmcglmm_global_pca_zscore_sensitivity.RData")

names(fullmcmc.results.global.pca.z) <- paste0("tree", 1:112)

epsnet.zscore.table.sens <-
    bind_rows(lapply(fullmcmc.results.global.pca.z,
                     function(x){res <- as.data.frame(x$Sol)
                         names(res) = c("Intercept", "Extinction Fraction", "NetDiversification", "PC1", "PC2", "NetDiversification:Extinction Fraction", "Extinction Fraction:PC1", "Extinction Fraction:PC2", "NetDiversification:PC1", "NetDiversification:PC2")
                         return(res)}))

epsnet.zscore.table.sens$tree <- rep(gsub("epsilon.", "", names(fulldata.net.bin)[26:137]), each = 1000)

epsnet.zscore.sens.plot <- melt(epsnet.zscore.table.sens)
## levels(epsnet.zscore.plot$variable) <- levels(epsnet.zscore.plot$variable)[c(1, 3, 2, 4:6, 9:10, 7:8)]

posterior.quants.epsnet.zscore.sens <- data.frame(
    variable = unique(epsnet.zscore.sens.plot$variable),
    aggregate(epsnet.zscore.sens.plot$value, by = list(epsnet.zscore.sens.plot$variable), FUN = function(x){quantile(x, c(0.025, 0.1, 0.5, 0.9, 0.975))})$x
)
names(posterior.quants.epsnet.zscore.sens) <- c("variable", "quant2.5", "quant10", "median", "quant90", "quant97.5")

medians.per.tree.epsnet.zscore.sens <- aggregate(epsnet.zscore.sens.plot$value, by = list(epsnet.zscore.sens.plot$tree, epsnet.zscore.sens.plot$variable), FUN = median)
names(medians.per.tree.epsnet.zscore.sens) <- c("tree", "variable", "value")

fulldata.net.bin$epsilon.z.sens <- unlist(aggregate(fulldata.net.bin$epsilon.eric0102, by = list(fulldata.net.bin$id), FUN = scale)$x)
fulldata.net.bin$netdiv.z.sens <- unlist(aggregate(fulldata.net.bin$netdiv.eric0102, by = list(fulldata.net.bin$id), FUN = scale)$x)
