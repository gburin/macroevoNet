library("tidyverse")
library("geosphere")
library("cowplot")
#source("../R/calculating_centrality.R")

elton.traits <- read.csv("../data/BirdFuncDat.txt", sep = "\t")

fulldata.net.bin$family <- elton.traits$BLFamilyLatin[match(fulldata.net.bin$species, gsub(" ", "_", elton.traits$Scientific))]
fulldata.net.bin$order <- elton.traits$IOCOrder[match(fulldata.net.bin$species, gsub(" ", "_", elton.traits$Scientific))]
fulldata.net.bin$genus <- sapply(as.character(fulldata.net.bin$species), function(x){strsplit(x, split = "_")[[1]][1]})

net.metadata <- net.metadata[-c(2, 7, 10, 18, 27),]

dist.nets <- distm(x = net.metadata[, c("long_dec", "lat_dec")], fun = distGeo)

mean.pca.family.net <- setNames(aggregate(fulldata.net.bin$pca, by = list(fulldata.net.bin$family, fulldata.net.bin$id), FUN = mean), c("family", "net", "pca"))
mean.pca.genus.net <- setNames(aggregate(fulldata.net.bin$pca, by = list(fulldata.net.bin$genus, fulldata.net.bin$id), FUN = mean), c("genus", "net", "pca"))

## Replacing NA's by 0 only in the network without the family

dist.correl.fam.single0 <- data.frame(dist = NA, correl = NA)    
for(i in 1:28){
    for(j in (i + 1):29){
        net1 <- mean.pca.family.net[mean.pca.family.net$net == unique(mean.pca.family.net$net)[i],]
        net2 <- mean.pca.family.net[mean.pca.family.net$net == unique(mean.pca.family.net$net)[j],]
        df <- data.frame(family = sort(unique(c(as.character(net1$family), as.character(net2$family)))),
                         net1 = NA,
                         net2 = NA)
        df$net1 <- net1$pca[match(as.character(df$family), net1$family)]
        df$net2 <- net2$pca[match(as.character(df$family), net2$family)]
        df$net1[is.na(df$net1)] <- 0
        df$net2[is.na(df$net2)] <- 0
        res <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                          correl = cor(df$net1, df$net2, method = "spearman")
                          )
        dist.correl.fam.single0 <- rbind(dist.correl.fam.single0, res)
    }
}
dist.correl.fam.single0 <- dist.correl.fam.single0[-1,]

dist.correl.gen.single0 <- data.frame(dist = NA, correl = NA)    
for(i in 1:28){
    for(j in (i + 1):29){
        net1 <- mean.pca.genus.net[mean.pca.genus.net$net == unique(mean.pca.genus.net$net)[i],]
        net2 <- mean.pca.genus.net[mean.pca.genus.net$net == unique(mean.pca.genus.net$net)[j],]
        df <- data.frame(genus = sort(unique(c(as.character(net1$genus), as.character(net2$genus)))),
                         net1 = NA,
                         net2 = NA)
        df$net1 <- net1$pca[match(as.character(df$genus), net1$genus)]
        df$net2 <- net2$pca[match(as.character(df$genus), net2$genus)]
        df$net1[is.na(df$net1)] <- 0
        df$net2[is.na(df$net2)] <- 0
        res <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                          correl = cor(df$net1, df$net2, method = "spearman")
                          )
        dist.correl.gen.single0 <- rbind(dist.correl.gen.single0, res)
    }
}
dist.correl.gen.single0 <- dist.correl.gen.single0[-1,]




## Replacing NA's by 0 in both networks when the family is absent in one of the networks

dist.correl.fam.full0 <- data.frame(dist = NA, correl = NA)    
for(i in 1:28){
    for(j in (i + 1):29){
        net1 <- mean.pca.family.net[mean.pca.family.net$net == unique(mean.pca.family.net$net)[i],]
        net2 <- mean.pca.family.net[mean.pca.family.net$net == unique(mean.pca.family.net$net)[j],]
        df <- data.frame(family = sort(unique(c(as.character(net1$family), as.character(net2$family)))),
                         net1 = NA,
                         net2 = NA)
        df$net1 <- net1$pca[match(as.character(df$family), net1$family)]
        df$net2 <- net2$pca[match(as.character(df$family), net2$family)]
        df[is.na(rowSums(df[, 2:3])), 2:3] <- c(0, 0)
        if(sum(df[, 2:3] != 0) > 0){
            res <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                              correl = cor(df$net1, df$net2, method = "spearman")
                              )
            dist.correl.fam.full0 <- rbind(dist.correl.fam.full0, res)
        } else {
            res <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                              correl = 0
                              )
            dist.correl.fam.full0 <- rbind(dist.correl.fam.full0, res)
        }
    }
}
dist.correl.fam.full0 <- dist.correl.fam.full0[-1,]

dist.correl.gen.full0 <- data.frame(dist = NA, correl = NA)    
for(i in 1:28){
    for(j in (i + 1):29){
        net1 <- mean.pca.genus.net[mean.pca.genus.net$net == unique(mean.pca.genus.net$net)[i],]
        net2 <- mean.pca.genus.net[mean.pca.genus.net$net == unique(mean.pca.genus.net$net)[j],]
        df <- data.frame(genus = sort(unique(c(as.character(net1$genus), as.character(net2$genus)))),
                         net1 = NA,
                         net2 = NA)
        df$net1 <- net1$pca[match(as.character(df$genus), net1$genus)]
        df$net2 <- net2$pca[match(as.character(df$genus), net2$genus)]
        df[is.na(rowSums(df[, 2:3])), 2:3] <- c(0, 0)
        if(sum(df[, 2:3] != 0) > 0){
            res <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                              correl = cor(df$net1, df$net2, method = "spearman")
                              )
            dist.correl.gen.full0 <- rbind(dist.correl.gen.full0, res)
        } else {
            res <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                              correl = 0
                              )
            dist.correl.gen.full0 <- rbind(dist.correl.gen.full0, res)
        }
    }
}
dist.correl.gen.full0 <- dist.correl.gen.full0[-1,]


load("../output/null_correl_geog_dist_Mar2021.RData")


family.within.single0 <-
    ggplot(data = dist.correl.fam.single0) +
    geom_line(data = null.smooth.family.within.single0, mapping = aes(x = x/1000, y = y, group = sim), colour = "lightgrey", alpha = 0.1) +
    geom_point(aes(x = dist/1000, y = correl), size = 1.5, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
    geom_smooth(aes(x = dist/1000, y = correl), se = FALSE, colour = "black", size = 2) +
    ylim(-1, 1) +
    labs(x = "Distance between nets (km)", y = "Correlation", colour = "No overlap", title = "Family - Within network") +
    theme(legend.position = "bottom") +
    theme_cowplot()

family.within.full0 <-
    ggplot(data = dist.correl.fam.full0) +
    geom_line(data = null.smooth.family.within.full0, mapping = aes(x = x/1000, y = y, group = sim), colour = "lightgrey", alpha = 0.1) +
    geom_point(aes(x = dist/1000, y = correl), size = 2, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
    geom_smooth(aes(x = dist/1000, y = correl), se = FALSE, colour = "black", size = 2) +
    ylim(-1, 1) +
    labs(x = "Distance between nets (km)", y = "Correlation", colour = "No overlap", title = "Family - Within network/full 0") +
    theme(legend.position = "bottom") +
    theme_cowplot()


genus.within.single0 <-
    ggplot(data = dist.correl.gen.single0) +
    geom_line(data = null.smooth.genus.within.single0, mapping = aes(x = x/1000, y = y, group = sim), colour = "lightgrey", alpha = 0.1) +
    geom_point(aes(x = dist/1000, y = correl), size = 1.5, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
    geom_smooth(aes(x = dist/1000, y = correl), se = FALSE, colour = "black", size = 2) +
    ylim(-1, 1) +
    labs(x = "Distance between nets (km)", y = "Correlation", colour = "No overlap", title = "Genus - Within network") +
    theme(legend.position = "bottom") +
    theme_cowplot()

genus.within.full0 <-
    ggplot(data = dist.correl.gen.full0) +
    geom_line(data = null.smooth.genus.within.full0, mapping = aes(x = x/1000, y = y, group = sim), colour = "lightgrey", alpha = 0.1) +
    geom_point(aes(x = dist/1000, y = correl), size = 2, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
    geom_smooth(aes(x = dist/1000, y = correl), se = FALSE, colour = "black", size = 2) +
    ylim(-1, 1) +
    labs(x = "Distance between nets (km)", y = "Correlation", colour = "No overlap", title = "Genus - Within network/full 0") +
    theme(legend.position = "bottom") +
    theme_cowplot()






family.between.single0 <-
    ggplot(data = dist.correl.fam.single0) +
    geom_line(data = null.smooth.family.between.single0, mapping = aes(x = x/1000, y = y, group = sim), colour = "lightgrey", alpha = 0.1) +
    geom_point(aes(x = dist/1000, y = correl), size = 1.5, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
    geom_smooth(aes(x = dist/1000, y = correl), se = FALSE, colour = "black", size = 2) +
    ylim(-1, 1) +
    labs(x = "Distance between nets (km)", y = "Correlation", colour = "No overlap", title = "Family - Between network") +
    theme(legend.position = "bottom") +
    theme_cowplot()

family.between.full0 <-
    ggplot(data = dist.correl.fam.full0) +
    geom_line(data = null.smooth.family.between.full0, mapping = aes(x = x/1000, y = y, group = sim), colour = "lightgrey", alpha = 0.1) +
    geom_point(aes(x = dist/1000, y = correl), size = 2, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
    geom_smooth(aes(x = dist/1000, y = correl), se = FALSE, colour = "black", size = 2) +
    ylim(-1, 1) +
    labs(x = "Distance between nets (km)", y = "Correlation", colour = "No overlap", title = "Family - Between network/full 0") +
    theme(legend.position = "bottom") +
    theme_cowplot()


genus.between.single0 <-
    ggplot(data = dist.correl.gen.single0) +
    geom_line(data = null.smooth.genus.between.single0, mapping = aes(x = x/1000, y = y, group = sim), colour = "lightgrey", alpha = 0.1) +
    geom_point(aes(x = dist/1000, y = correl), size = 1.5, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
    geom_smooth(aes(x = dist/1000, y = correl), se = FALSE, colour = "black", size = 2) +
    ylim(-1, 1) +
    labs(x = "Distance between nets (km)", y = "Correlation", colour = "No overlap", title = "Genus - Between network") +
    theme(legend.position = "bottom") +
    theme_cowplot()

genus.between.full0 <-
    ggplot(data = dist.correl.gen.full0) +
    geom_line(data = null.smooth.genus.between.full0, mapping = aes(x = x/1000, y = y, group = sim), colour = "lightgrey", alpha = 0.1) +
    geom_point(aes(x = dist/1000, y = correl), size = 2, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
    geom_smooth(aes(x = dist/1000, y = correl), se = FALSE, colour = "black", size = 2) +
    ylim(-1, 1) +
    labs(x = "Distance between nets (km)", y = "Correlation", colour = "No overlap", title = "Genus - Between network/full 0") +
    theme(legend.position = "bottom") +
    theme_cowplot()




plot_grid(family.within.single0,
          family.within.full0,
          genus.within.single0,
          genus.within.full0,
          nrow = 2,
          align = 'hv',
          labels = LETTERS[1:4]
          )


plot_grid(family.between.single0,
          family.between.full0,
          genus.between.single0,
          genus.between.full0,
          nrow = 2,
          align = 'hv',
          labels = LETTERS[1:4]
          )

fig3 <-
    plot_grid(
        family.within.single0,
        family.between.single0,
        genus.within.single0,
        genus.between.single0,
        nrow = 2,
        align = 'hv',
        labels = LETTERS[1:4]
    )

#ggsave("../output/fig3_final.pdf", plot = fig3, width = 11, height = 7, units = "in")
