library("MCMCglmm")
library("ape")
library("tidyverse")
library("reshape2")
library("cowplot")
library("RColorBrewer")

### Global + PCA

load("../output/mcmcglmm_global_pca_zscore_sensitivity_lambdamu.RData")

names(fullmcmc.results.global.pca.z) <- paste0("tree", 1:112)

lambdamu.zscore.table.sens <-
    bind_rows(lapply(fullmcmc.results.global.pca.z,
                     function(x){res <- as.data.frame(x$Sol); names(res) = c("Intercept", "Speciation", "Extinction", "PC1", "PC2", "Speciation:Extinction", "Speciation:PC1", "Speciation:PC2", "Extinction:PC1", "Extinction:PC2")
                         return(res)}))

lambdamu.zscore.table.sens$tree <- rep(gsub("epsilon.", "", names(fulldata.net.bin)[26:137]), each = 1000)

lambdamu.zscore.sens.plot <- melt(lambdamu.zscore.table.sens)
levels(lambdamu.zscore.sens.plot$variable) <- levels(lambdamu.zscore.sens.plot$variable)[c(1:6, 7:10)]

posterior.quants.lambdamu.zscore.sens <- data.frame(
    variable = unique(lambdamu.zscore.sens.plot$variable),
    aggregate(lambdamu.zscore.sens.plot$value, by = list(lambdamu.zscore.sens.plot$variable), FUN = function(x){quantile(x, c(0.025, 0.1, 0.5, 0.9, 0.975))})$x
)
names(posterior.quants.lambdamu.zscore.sens) <- c("variable", "quant2.5", "quant10", "median", "quant90", "quant97.5")

medians.per.tree.lambdamu.zscore.sens <- aggregate(lambdamu.zscore.sens.plot$value, by = list(lambdamu.zscore.sens.plot$tree, lambdamu.zscore.sens.plot$variable), FUN = median)
names(medians.per.tree.lambdamu.zscore.sens) <- c("tree", "variable", "value")

fulldata.net.bin$lambda.z.sens <- unlist(aggregate(fulldata.net.bin$lambda.eric0102, by = list(fulldata.net.bin$id), FUN = scale)$x)
fulldata.net.bin$mu.z.sens <- unlist(aggregate(fulldata.net.bin$mu.eric0102, by = list(fulldata.net.bin$id), FUN = scale)$x)
