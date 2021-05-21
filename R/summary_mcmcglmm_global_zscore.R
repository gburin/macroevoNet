library("MCMCglmm")
library("ape")
library("tidyverse")
library("reshape2")
library("cowplot")
library("RColorBrewer")

### Global + PCA

load("../output/mcmcglmm_global_pca_zscore.RData")

names(fullmcmc.results.global.pca.z) <- paste0("tree", 1:112)

epsnet.zscore.table <-
    bind_rows(lapply(fullmcmc.results.global.pca.z,
                     function(x){res <- as.data.frame(x$Sol)
                         names(res) = c("Intercept", "Epsilon", "NetDiversification", "PC1", "PC2", "NetDiversification:Epsilon", "Epsilon:PC1", "Epsilon:PC2:PC2", "NetDiversification:PC1", "NetDiversification:PC2")
                         return(res)}))

epsnet.zscore.table$tree <- rep(gsub("epsilon.", "", names(fulldata.net.bin)[26:137]), each = 1000)

epsnet.zscore.plot <- melt(epsnet.zscore.table)
## levels(epsnet.zscore.plot$variable) <- levels(epsnet.zscore.plot$variable)[c(1, 3, 2, 4:6, 9:10, 7:8)]

posterior.quants.epsnet.zscore <- data.frame(
    variable = unique(epsnet.zscore.plot$variable),
    aggregate(epsnet.zscore.plot$value, by = list(epsnet.zscore.plot$variable), FUN = function(x){quantile(x, c(0.025, 0.1, 0.5, 0.9, 0.975))})$x
)
names(posterior.quants.epsnet.zscore) <- c("variable", "quant2.5", "quant10", "median", "quant90", "quant97.5")

medians.per.tree.epsnet.zscore <- aggregate(epsnet.zscore.plot$value, by = list(epsnet.zscore.plot$tree, epsnet.zscore.plot$variable), FUN = median)
names(medians.per.tree.epsnet.zscore) <- c("tree", "variable", "value")

fulldata.net.bin$epsilon.z <- unlist(aggregate(fulldata.net.bin$epsilon.eric0102, by = list(fulldata.net.bin$id), FUN = scale)$x)
fulldata.net.bin$netdiv.z <- unlist(aggregate(fulldata.net.bin$netdiv.eric0102, by = list(fulldata.net.bin$id), FUN = scale)$x)
