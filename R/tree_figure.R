library("ggtree")
library("ape")
library("ggsci")
library("ape")
library("phytools")
library("ggnewscale")
library("viridis")
library("RColorBrewer")
library("ggforce")
library("ggstance")
library("reshape2")

source("../R/calculating_centrality.R")

## Importing Rates
rates.eric <- read.csv("../data/avg_rates_eric_bamm.csv")
rates.hack <- read.csv("../data/avg_rates_hack_bamm.csv")

rates <- rbind(rates.eric, rates.hack)

## Reorganizing tables
lambda <- rates[rates$param == "lambda",]
mu <- rates[rates$param == "mu",]

epsilon <- mu
epsilon$value  <- epsilon$value / lambda$value
netdiv <- lambda
netdiv$value <- netdiv$value - mu$value

lambda.split <- reshape2::dcast(lambda[,-4], species ~ .id)
mu.split <- reshape2::dcast(mu[,-4], species ~ .id)
epsilon.split <- reshape2::dcast(epsilon[,-4], species ~ .id)
netdiv.split <- reshape2::dcast(netdiv[,-4], species ~ .id)

fulldata.net.bin <- cbind(fulldata.net.bin, setNames(epsilon.split[match(fulldata.net.bin$species, epsilon.split$species), -1], paste0("epsilon.", names(epsilon.split)[-1])))
fulldata.net.bin <- cbind(fulldata.net.bin, setNames(netdiv.split[match(fulldata.net.bin$species, netdiv.split$species), -1], paste0("netdiv.", names(netdiv.split)[-1])))
fulldata.net.bin <- cbind(fulldata.net.bin, setNames(lambda.split[match(fulldata.net.bin$species, lambda.split$species), -1], paste0("lambda.", names(lambda.split)[-1])))
fulldata.net.bin <- cbind(fulldata.net.bin, setNames(mu.split[match(fulldata.net.bin$species, mu.split$species), -1], paste0("mu.", names(mu.split)[-1])))

## Removing species without rates

fulldata.net.bin <- fulldata.net.bin[!is.na(fulldata.net.bin$epsilon.eric0102),]

clim.pca <- prcomp(net.metadata[, 14:17], scale = TRUE)

fulldata.net.bin <- cbind(fulldata.net.bin, clim.pca$x[match(fulldata.net.bin$id, net.metadata$Network),])
names(fulldata.net.bin)[(ncol(fulldata.net.bin) - 3):ncol(fulldata.net.bin)] <- paste0("clim.pc", 1:4)

tree <- read.tree("../data/full_bird_tree.txt")
elton.traits <- read.csv("../data/BirdFuncDat.txt", sep = "\t")

fulldata.net.bin$family <- as.character(elton.traits$BLFamilyLatin[match(fulldata.net.bin$species, gsub(" ", "_", elton.traits$Scientific))])
fulldata.net.bin$order <- as.character(elton.traits$IOCOrder[match(fulldata.net.bin$species, gsub(" ", "_", elton.traits$Scientific))])

tr.drop <- drop.tip(tree, tip = tree$tip.label[is.na(match(tree$tip.label, fulldata.net.bin$species))])



order.nodes <- sapply(unique(fulldata.net.bin$order), function(x){ape::getMRCA(ladderize(tr.drop), as.character(fulldata.net.bin$species[fulldata.net.bin$order == x]))})

family.nodes <- sapply(unique(fulldata.net.bin$family), function(x){ape::getMRCA(ladderize(tr.drop), as.character(fulldata.net.bin$species[fulldata.net.bin$family == x]))})
family.nodes <- family.nodes[-sapply(family.nodes, is.null) == FALSE]
clade.cols <- pal_npg(alpha = 1)
family.cols <- colorRampPalette(pal_lancet(alpha = 1)(9))

data.plot <- fulldata.net.bin[, c(2, 9, 26, 138, 250, 362)]
names(data.plot)[1] <- "id"

data.plot$nnet <- as.integer(table(as.character(data.plot$id))[match(data.plot$id, names(table(as.character(data.plot$id))))])



tree.bars <-
    ggtree(ladderize(tr.drop)) +
    geom_cladelabel(node = order.nodes[1], color = clade.cols(9)[1], label = substr(names(order.nodes[1]), 1, 3), barsize = 3, angle = 45) +
    geom_cladelabel(node = order.nodes[2], color = clade.cols(9)[2], label = substr(names(order.nodes[2]), 1, 3), barsize = 3, angle = 45) +
    geom_cladelabel(node = order.nodes[3], color = clade.cols(9)[3], label = substr(names(order.nodes[3]), 1, 3), barsize = 3, angle = 45) +
    geom_cladelabel(node = order.nodes[4], color = clade.cols(9)[4], label = substr(names(order.nodes[4]), 1, 3), barsize = 3, angle = 45) +
    geom_cladelabel(node = order.nodes[5], color = clade.cols(9)[5], label = substr(names(order.nodes[5]), 1, 3), barsize = 3, angle = 45) +
    geom_cladelabel(node = order.nodes[6], color = clade.cols(9)[6], label = substr(names(order.nodes[6]), 1, 3), barsize = 3, angle = 45) +
    geom_cladelabel(node = order.nodes[7], color = clade.cols(9)[7], label = substr(names(order.nodes[7]), 1, 3), barsize = 3, angle = 45) +
    geom_cladelabel(node = order.nodes[8], color = clade.cols(9)[8], label = substr(names(order.nodes[8]), 1, 3), barsize = 3, angle = 45) +
    geom_cladelabel(node = order.nodes[9], color = clade.cols(9)[9], label = substr(names(order.nodes[9]), 1, 3), barsize = 3, angle = 45)


p1 <-
    facet_plot(tree.bars, panel = "Speciation", data = data.plot, geom = geom_tile, mapping = aes(x = 1, fill = lambda.eric0102, group = label), colour = "white", width = 0.1) +
    scale_fill_gradient(low = paste0(brewer.pal(4, "Dark2")[2], "25"), high = paste0(brewer.pal(4, "Dark2")[2], "FF"), name = "Diversification\nRate") +
    theme(legend.position = "none")

p2 <- p1 + new_scale_fill()

p3 <-
    facet_plot(p2, panel = "Extinction", data = data.plot, geom = geom_tile, mapping = aes(x = 1, fill = mu.eric0102, group = label), colour = "white", width = 0.1) +
    scale_fill_gradient(low = paste0(brewer.pal(4, "Dark2")[4], "25"), high = paste0(brewer.pal(4, "Dark2")[4], "FF"), name = "Diversification\nRate") +
    theme(legend.position = "none")

p4 <- p3 + new_scale_fill()

p5 <-
    facet_plot(p4, panel = "Extinction Fraction", data = data.plot, geom = geom_tile, mapping = aes(x = 1, fill = epsilon.eric0102, group = label), colour = "white", width = 0.1) +
    scale_fill_gradient(low = paste0(brewer.pal(3, "Dark2")[3], "25"), high = paste0(brewer.pal(3, "Dark2")[3], "FF"), name = "Diversification\nRate") +
    theme(legend.position = "none")

p6 <- p5 + new_scale_fill()

p7 <-
    facet_plot(p6, panel = "Net Diversification", data = data.plot, geom = geom_tile, mapping = aes(x = 1, fill = netdiv.eric0102, group = label), colour = "white", width = 0.1) +
    scale_fill_gradient(low = paste0(brewer.pal(3, "Dark2")[1], "25"), high = paste0(brewer.pal(3, "Dark2")[1], "FF"), name = "Diversification\nRate") +
    theme(legend.position = "none")

#p4 <- p3 + new_scale_color()

p8 <-
    facet_plot(p7, panel = "Centrality", data = data.plot, geom = geom_boxploth, mapping = aes(x = pca, group = label, colour = factor(nnet)), outlier.size = 0.5) +
    scale_colour_manual(values = rev(brewer.pal(6, "Spectral"))) +
    theme_tree2()

#facet_widths(p5, c(0.6, 0.05, 0.05, 0.3))

ggsave("../output/tree_fig_full.pdf", p8)

