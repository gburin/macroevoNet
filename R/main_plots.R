library("rmarkdown")
library("ggfortify")
library("ggdist")
library("cowplot")
library("ggsci")
library("patchwork")
source("../R/data_preparation_mcmcglmm.R")
source("../R/summary_mcmcglmm_global_zscore.R")
source("../R/summary_mcmcglmm_lambdamu_zscore.R")
source("../R/geog_similarity_plots.R")
source("../R/plots_diet.R")
source("../R/tree_figure.R")


## Figure 2

## Workaround for plotting pointinterval

epsnet.zscore.plot$pi.x <- NA
epsnet.zscore.plot$pi.x <- seq(0, 10)[match(epsnet.zscore.plot$variable, levels(epsnet.zscore.plot$variable)[c(1, 3, 2, 4:6, 9:10, 7:8)])]

lambdamu.zscore.plot$pi.x <- NA
lambdamu.zscore.plot$pi.x <- seq(0, 10)[match(lambdamu.zscore.plot$variable, levels(lambdamu.zscore.plot$variable))]

fulldata.net.bin$lambda.z <- unlist(aggregate(fulldata.net.bin$lambda.eric0102, by = list(fulldata.net.bin$id), FUN = scale)$x)
fulldata.net.bin$mu.z <- unlist(aggregate(fulldata.net.bin$mu.eric0102, by = list(fulldata.net.bin$id), FUN = scale)$x)
fulldata.net.bin$epsilon.z <- unlist(aggregate(fulldata.net.bin$epsilon.eric0102, by = list(fulldata.net.bin$id), FUN = scale)$x)
fulldata.net.bin$netdiv.z <- unlist(aggregate(fulldata.net.bin$netdiv.eric0102, by = list(fulldata.net.bin$id), FUN = scale)$x)

## post.main <-
##     ggplot() +
##     geom_hline(yintercept = 0, linetype = "dashed", colour = "lightgrey") +
##     stat_slab(data = subset(epsnet.zscore.plot, variable != "Intercept"), aes(x = factor(variable, levels = levels(variable)[c(1, 3, 2, 4:6, 9:10, 7:8)]), y = value, colour = factor(variable, levels = levels(variable)[c(1, 3, 2, 4:6, 9:10, 7:8)])), position = "dodge", fill = NA, alpha = 0.2) +
##     stat_slab(data = subset(lambdamu.zscore.plot, variable != "Intercept"), aes(x = variable, y = value, colour = variable, fill = variable), side = "left", position = "dodge", alpha = 0.2) +
##     stat_pointinterval(data = subset(epsnet.zscore.plot, variable != "Intercept"), aes(x = 1 + 0.5, y = value, colour = factor(variable, levels = levels(variable)[c(1, 3, 2, 4:6, 9:10, 7:8)])), position = "dodge", interval_size_range = c(1.5, 3)) +
##     stat_pointinterval(data = subset(lambdamu.zscore.plot, variable != "Intercept"), aes(x = 1 - 0.5, y = value, colour = variable), position = "dodge", interval_size_range = c(1.5, 3)) +
##     geom_point(data = subset(medians.per.tree.epsnet.zscore, variable != "Intercept"), aes(x = 0.95, y = value, colour = factor(variable, levels = levels(variable)[c(1, 3, 2, 4:6, 9:10, 7:8)])), shape = 19, size = 0.75, alpha = 0.3) +
##     geom_point(data = subset(medians.per.tree.lambdamu.zscore, variable != "Intercept"), aes(x = 1.05, y = value, colour = variable), shape = 19, size = 0.75, alpha = 0.3) +
##     #labs(title = "stat_gradientinterval(position = 'dodge')") +
##     scale_colour_aaas(alpha = 0.7) +
##     scale_fill_aaas(alpha = 0.7) +
##     labs(y = "Posterior Density") +
##     theme_cowplot() +
##     theme(legend.position = "none",
##           axis.title.x = element_blank(),
##           axis.text.x = element_blank(),
##           axis.ticks.x = element_blank(),
##           axis.line.x = element_blank(),
##           strip.placement = "inside",
##           strip.background = element_blank()) +
##     facet_wrap(. ~ factor(variable, levels = levels(variable)[c(1, 3, 2, 4:6, 9:10, 7:8)]), ncol = 9, scales = "free_x", strip.position = "bottom")

epsnetdiv.grad <-
    ggplot() +
    geom_point(data = subset(medians.per.tree.epsnet.zscore, variable != "Intercept"), aes(x = factor(variable, levels = levels(variable)[c(1, 3, 2, 4:6, 9:10, 7:8)]), y = value, colour = factor(variable, levels = levels(variable)[c(1, 3, 2, 4:6, 9:10, 7:8)])), shape = 19, size = 0.5) +
    stat_gradientinterval(data = subset(epsnet.zscore.plot, variable != "Intercept"), aes(x = factor(variable, levels = levels(variable)[c(1, 3, 2, 4:6, 9:10, 7:8)]), y = value, colour = factor(variable, levels = levels(variable)[c(1, 3, 2, 4:6, 9:10, 7:8)]), fill = factor(variable, levels = levels(variable)[c(1, 3, 2, 4:6, 9:10, 7:8)]))) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "lightgrey") +
    scale_colour_aaas(alpha = 1) +
    scale_fill_aaas(alpha = 0.01) +
    ylim(-1.5, 1.5) +
    labs(y = "Posterior Density") +
    theme_cowplot() +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())

lambdamu.grad <-
    ggplot() +
    geom_point(data = subset(medians.per.tree.lambdamu.zscore, variable != "Intercept"), aes(x = variable, y = value, colour = variable), shape = 19, size = 0.5) +
    stat_gradientinterval(data = subset(lambdamu.zscore.plot, variable != "Intercept"), aes(x = variable, y = value, colour = variable, fill = variable)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "lightgrey") +
    scale_colour_aaas(alpha = 1) +
    scale_fill_aaas(alpha = 0.01) +
    ylim(-1.5, 1.5) +
    labs(y = "Posterior Density") +
    theme_cowplot() +
    theme(legend.position = "none", axis.title.x = element_blank())

ribbon.lambda <- tibble(
    intercept = lambdamu.zscore.table[, 1],
    slope = lambdamu.zscore.table[, 2],
    x = list(seq(min(fulldata.net.bin$lambda.z), max(fulldata.net.bin$lambda.z), by = 0.1)),
    y = map2(intercept, slope, ~ .x + .y * seq(min(fulldata.net.bin$lambda.z), max(fulldata.net.bin$lambda.z), by = 0.1))
) %>%
    unnest(c(x, y))

ribbon.lambda.plot <-
    ribbon.lambda %>%
    group_by(x) %>%
    median_qi(y, .width = c(0.5, 0.8, 0.95))

ribbon.mu <- tibble(
    intercept = lambdamu.zscore.table[, 1],
    slope = lambdamu.zscore.table[, 3],
    x = list(seq(min(fulldata.net.bin$mu.z), max(fulldata.net.bin$mu.z), by = 0.1)),
    y = map2(intercept, slope, ~ .x + .y * seq(min(fulldata.net.bin$mu.z), max(fulldata.net.bin$mu.z), by = 0.1))
) %>%
    unnest(c(x, y))

ribbon.mu.plot <-
    ribbon.mu %>%
    group_by(x) %>%
    median_qi(y, .width = c(0.5, 0.8, 0.95))

ribbon.epsilon <- tibble(
    intercept = epsnet.zscore.table[, 1],
    slope = epsnet.zscore.table[, 2],
    x = list(seq(min(fulldata.net.bin$epsilon.z), max(fulldata.net.bin$epsilon.z), by = 0.1)),
    y = map2(intercept, slope, ~ .x + .y * seq(min(fulldata.net.bin$epsilon.z), max(fulldata.net.bin$epsilon.z), by = 0.1))
) %>%
    unnest(c(x, y))

ribbon.epsilon.plot <-
    ribbon.epsilon %>%
    group_by(x) %>%
    median_qi(y, .width = c(0.5, 0.8, 0.95))

ribbon.netdiv <- tibble(
    intercept = epsnet.zscore.table[, 1],
    slope = epsnet.zscore.table[, 3],
    x = list(seq(min(fulldata.net.bin$netdiv.z), max(fulldata.net.bin$netdiv.z), by = 0.1)),
    y = map2(intercept, slope, ~ .x + .y * seq(min(fulldata.net.bin$netdiv.z), max(fulldata.net.bin$netdiv.z), by = 0.1))
) %>%
    unnest(c(x, y))

ribbon.netdiv.plot <-
    ribbon.netdiv %>%
    group_by(x) %>%
    median_qi(y, .width = c(0.5, 0.8, 0.95))


lambda.main <-
    ggplot(data = ribbon.lambda.plot) +
    geom_lineribbon(aes(x = x, y = y, ymin = .lower, ymax = .upper, color = pal_aaas(alpha = 1)(2)[1])) +
    geom_segment(aes(x = min(x), xend = max(x), y = min(y), yend = max(y)), color = pal_aaas(alpha = 1)(2)[1], size = 1.5) +
    geom_point(data = fulldata.net.bin, aes(x = lambda.z, y = pca), alpha = 0.3) +
    scale_fill_manual(values = c(pal_aaas(alpha = 0.1)(2)[1], pal_aaas(alpha = 0.1)(2)[1], pal_aaas(alpha = 0.1)(2)[1])) +
    labs(x = "Speciation (Standardized)", y = "PCA centrality") +
    ylim(-4, 10) +
    theme_cowplot() +
    theme(legend.position = "none")

mu.main <-
    ggplot(data = ribbon.mu.plot) +
    geom_lineribbon(aes(x = x, y = y, ymin = .lower, ymax = .upper)) +
    geom_segment(aes(x = min(x), xend = max(x), y = max(y), yend = min(y)), color = pal_aaas(alpha = 1)(2)[2], size = 1.5) +
    geom_point(data = fulldata.net.bin, aes(x = mu.z, y = pca), alpha = 0.3) +
    scale_fill_manual(values = c(pal_aaas(alpha = 0.1)(2)[2], pal_aaas(alpha = 0.1)(2)[2], pal_aaas(alpha = 0.1)(2)[2])) +
    labs(x = "Extinction (Standardized)", y = element_blank()) +
    ylim(-4, 10) +
    theme_cowplot() +
    theme(legend.position = "none")

epsilon.main <-
    ggplot(data = ribbon.epsilon.plot) +
    geom_lineribbon(aes(x = x, y = y, ymin = .lower, ymax = .upper)) +
    geom_segment(aes(x = min(x), xend = max(x), y = max(y), yend = min(y)), color = pal_aaas(alpha = 1)(2)[2], size = 1.5) +
    geom_point(data = fulldata.net.bin, aes(x = epsilon.z, y = pca), alpha = 0.3) +
    scale_fill_manual(values = c(pal_aaas(alpha = 0.1)(2)[2], pal_aaas(alpha = 0.1)(2)[2], pal_aaas(alpha = 0.1)(2)[2])) +
    labs(x = "Extinction Fraction (Standardized)", y = element_blank()) +
    ylim(-4, 10) +
    theme_cowplot() +
    theme(legend.position = "none")

netdiv.main <-
    ggplot(data = ribbon.netdiv.plot) +
    geom_lineribbon(aes(x = x, y = y, ymin = .lower, ymax = .upper)) +
    geom_segment(aes(x = min(x), xend = max(x), y = min(y), yend = max(y)), color = pal_aaas(alpha = 1)(2)[1], size = 1.5) +
    geom_point(data = fulldata.net.bin, aes(x = netdiv.z, y = pca), alpha = 0.3) +
    scale_fill_manual(values = c(pal_aaas(alpha = 0.1)(2)[1], pal_aaas(alpha = 0.1)(2)[1], pal_aaas(alpha = 0.1)(2)[1])) +
    labs(x = "Net Diversification (Standardized)", y = element_blank()) +
    ylim(-4, 10) +
    theme_cowplot() +
    theme(legend.position = "none")

fig2.main <- ((lambda.main | mu.main | netdiv.main| epsilon.main ) / post.main) +
    plot_layout(heights = c(0.4, 1))

ggsave("../output/fig2_full.svg", plot = fig2.main, dpi = 600)

fig2.alternative <-
    ((lambda.main | mu.main | netdiv.main| epsilon.main ) / (lambdamu.grad + epsnetdiv.grad)) +
    plot_layout(heights = c(0.4, 1))

ggsave("../output/temp_figs/fig2_full_alternative.svg", plot = fig2.alternative, dpi = 600, width = 20, height = 10)

## Figure 3

int.tax <- read.csv("../data/interaction_per_taxonomy.csv")

pc.diet <-
    ggplot(data = pcs.diet.mass) +
    geom_point(aes(x = PC1, y = PC2, colour = diet), size = 3, alpha = 0.3) +
    scale_colour_d3() +
    theme(legend.position = "bottom") +
    geom_segment(x = 0, xend = pPCA.diet.mass$Evec[1,1], y = 0, yend = pPCA.diet.mass$Evec[1,2], colour = brewer.pal(3, "Set1")[1], arrow = arrow(length = unit(0.03, "npc"))) +
    geom_text(x = pPCA.diet.mass$Evec[1,1] + 0.05, y = pPCA.diet.mass$Evec[1,2] - 0.05, label = "Invertebrates", colour = brewer.pal(3, "Set1")[1], size = 7) +
    geom_segment(x = 0, xend = pPCA.diet.mass$Evec[5,1], y = 0, yend = pPCA.diet.mass$Evec[5,2], colour = brewer.pal(3, "Set1")[1], arrow = arrow(length = unit(0.03, "npc"))) +
    geom_text(x = pPCA.diet.mass$Evec[5,1] - 0.05, y = pPCA.diet.mass$Evec[5,2] - 0.05, label = "Fruits", colour = brewer.pal(3, "Set1")[1], size = 7) +
    geom_segment(x = 0, xend = pPCA.diet.mass$Evec[7,1], y = 0, yend = pPCA.diet.mass$Evec[7,2], colour = brewer.pal(3, "Set1")[1], arrow = arrow(length = unit(0.03, "npc"))) +
    geom_text(x = pPCA.diet.mass$Evec[7,1], y = pPCA.diet.mass$Evec[7,2] + 0.05, label = "Seeds", colour = brewer.pal(3, "Set1")[1], size = 7) +
    labs(x = "PC1 (51.4%)", y = "PC2 (28.5%", colour = "Dietary Guild") +
    ylim(-0.5, 1.2) +
    xlim(-1, 1) +
    theme_cowplot() +
    theme(legend.position = "top")

pccent.phylodiv <-
    ggplot(data = int.tax) +
    geom_point(aes(x = pca, y = phylodiv)) +
    geom_smooth(aes(x = pca, y = phylodiv), method = "lm", se = FALSE, colour = "red") +
    labs(x = "PCA centrality", y = "Phylogenetic Diversity of Plants") +
    theme_cowplot()

upper.row <-
    plot_grid(
        pccent.phylodiv,
        pc.diet,
        ncol = 2,
        align = 'v',
        labels = LETTERS[1:2]
    )

lower.row <-
    plot_grid(disp.1.av.pc1,
              disp.1.mnnd.pc1,
              ncol = 2,
              align = 'hv',
              labels = LETTERS[3:4]
              )

## fig4 <- plot_grid(upper.row,
##           lower.row,
##           nrow = 2,
##     #      rel_widths = c(0.6, 0.4),
##           align = 'h'
##           )

fig4 <-
    plot_grid(
        disp.1.av.pc1,
        disp.1.mnnd.pc1,
        pccent.phylodiv,
        ncol = 3,
        align = 'h',
        labels = LETTERS[1:3]
    )


## Figure 4

plot_grid(family.within.single0,
          genus.within.single0,
          family.within.full0,
          genus.within.full0,
          nrow = 2,
          align = 'hv',
          labels = LETTERS[1:4]
          )
