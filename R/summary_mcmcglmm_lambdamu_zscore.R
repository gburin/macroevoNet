library("MCMCglmm")
library("ape")
library("tidyverse")
library("reshape2")
library("cowplot")
library("RColorBrewer")

### Global + PCA

load("../output/mcmcglmm_lambdamu_pca_zscore.RData")

names(fullmcmc.results.lambdamu.pca.z1) <- paste0("tree", 1:112)

lambdamu.zscore.table <-
    bind_rows(lapply(fullmcmc.results.lambdamu.pca.z1,
                     function(x){res <- as.data.frame(x$Sol); names(res) = c("Intercept", "Speciation", "Extinction", "PC1", "PC2", "Speciation:Extinction", "Speciation:PC1", "Speciation:PC2", "Extinction:PC1", "Extinction:PC2")
                         return(res)}))

lambdamu.zscore.table$tree <- rep(gsub("epsilon.", "", names(fulldata.net.bin)[26:137]), each = 1000)

lambdamu.zscore.plot <- melt(lambdamu.zscore.table)
levels(lambdamu.zscore.plot$variable) <- levels(lambdamu.zscore.plot$variable)[c(1:6, 7:10)]

posterior.quants.lambdamu.zscore <- data.frame(
    variable = unique(lambdamu.zscore.plot$variable),
    aggregate(lambdamu.zscore.plot$value, by = list(lambdamu.zscore.plot$variable), FUN = function(x){quantile(x, c(0.025, 0.1, 0.5, 0.9, 0.975))})$x
)
names(posterior.quants.lambdamu.zscore) <- c("variable", "quant2.5", "quant10", "median", "quant90", "quant97.5")

medians.per.tree.lambdamu.zscore <- aggregate(lambdamu.zscore.plot$value, by = list(lambdamu.zscore.plot$tree, lambdamu.zscore.plot$variable), FUN = median)
names(medians.per.tree.lambdamu.zscore) <- c("tree", "variable", "value")

fulldata.net.bin$lambda.z <- unlist(aggregate(fulldata.net.bin$lambda.eric0102, by = list(fulldata.net.bin$id), FUN = scale)$x)
fulldata.net.bin$mu.z <- unlist(aggregate(fulldata.net.bin$mu.eric0102, by = list(fulldata.net.bin$id), FUN = scale)$x)
