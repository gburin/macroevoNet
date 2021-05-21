library("tidyverse")
library("geosphere")
library("foreach")
library("doMC")
registerDoMC(56)
source("../R/calculating_centrality.R")

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





### Randomizing within networks

rand.geog.nets.single0 <- function(x, level = "family"){
    res <- data.frame(dist = NA, correl = NA)
    for(i in 1:28){
        for(j in (i + 1):29){
            net1.full <- fulldata.net.bin[fulldata.net.bin$id == unique(fulldata.net.bin$id)[i],]
            net2.full <- fulldata.net.bin[fulldata.net.bin$id == unique(fulldata.net.bin$id)[j],]
            net1.full$pca <- sample(net1.full$pca)
            net2.full$pca <- sample(net2.full$pca)
            if(level == "family"){
                mean.pca.family.net1 <- setNames(aggregate(net1.full$pca, by = list(net1.full$family), FUN = mean), c("family", "pca"))
                mean.pca.family.net2 <- setNames(aggregate(net2.full$pca, by = list(net2.full$family), FUN = mean), c("family", "pca"))
                df <- data.frame(family = sort(unique(c(as.character(mean.pca.family.net1$family), as.character(mean.pca.family.net2$family)))),
                                 net1 = NA,
                                 net2 = NA)
                df$net1 <- mean.pca.family.net1$pca[match(as.character(df$family), as.character(mean.pca.family.net1$family))]
                df$net2 <- mean.pca.family.net2$pca[match(as.character(df$family), as.character(mean.pca.family.net2$family))]
                df$net1[is.na(df$net1)] <- 0
                df$net2[is.na(df$net2)] <- 0
                res.df <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                                     correl = cor(df$net1, df$net2, method = "spearman")
                                     )
                res <- rbind(res, res.df)
            } else {
                net1.full$genus <- sapply(net1.full$species, function(x){strsplit(as.character(x), split = "_")[[1]][1]})
                net2.full$genus <- sapply(net2.full$species, function(x){strsplit(as.character(x), split = "_")[[1]][1]})
                mean.pca.genus.net1 <- setNames(aggregate(net1.full$pca, by = list(net1.full$genus), FUN = mean), c("genus", "pca"))
                mean.pca.genus.net2 <- setNames(aggregate(net2.full$pca, by = list(net2.full$genus), FUN = mean), c("genus", "pca"))
                df <- data.frame(genus = sort(unique(c(as.character(mean.pca.genus.net1$genus), as.character(mean.pca.genus.net2$genus)))),
                                 net1 = NA,
                                 net2 = NA)
                df$net1 <- mean.pca.genus.net1$pca[match(as.character(df$genus), as.character(mean.pca.genus.net1$genus))]
                df$net2 <- mean.pca.genus.net2$pca[match(as.character(df$genus), as.character(mean.pca.genus.net2$genus))]
                df$net1[is.na(df$net1)] <- 0
                df$net2[is.na(df$net2)] <- 0
                res.df <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                                     correl = cor(df$net1, df$net2, method = "spearman")
                                     )
                res <- rbind(res, res.df)
            }
        }
    }
    res <- res[-1,]
    rand.dist.correl <-
        ggplot(data = res) +
        geom_point(aes(x = dist, y = correl), size = 2) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
        geom_smooth(aes(x = dist, y = correl), se = FALSE) +
        ylim(-1, 1)
    rand.dist.correl.build <- ggplot_build(rand.dist.correl)$data[[3]][, c("x", "y")]
    return(list(sims = res, data.plot = rand.dist.correl.build))
}


rand.geog.nets.full0 <- function(x, level = "family"){
    res <- data.frame(dist = NA, correl = NA)
    for(i in 1:28){
        for(j in (i + 1):29){
            net1.full <- fulldata.net.bin[fulldata.net.bin$id == unique(fulldata.net.bin$id)[i],]
            net2.full <- fulldata.net.bin[fulldata.net.bin$id == unique(fulldata.net.bin$id)[j],]
            net1.full$pca <- sample(net1.full$pca)
            net2.full$pca <- sample(net2.full$pca)
            if(level == "family"){
                mean.pca.family.net1 <- setNames(aggregate(net1.full$pca, by = list(net1.full$family), FUN = mean), c("family", "pca"))
                mean.pca.family.net2 <- setNames(aggregate(net2.full$pca, by = list(net2.full$family), FUN = mean), c("family", "pca"))
                df <- data.frame(family = sort(unique(c(as.character(mean.pca.family.net1$family), as.character(mean.pca.family.net2$family)))),
                                 net1 = NA,
                                 net2 = NA)
                df$net1 <- mean.pca.family.net1$pca[match(as.character(df$family), as.character(mean.pca.family.net1$family))]
                df$net2 <- mean.pca.family.net2$pca[match(as.character(df$family), as.character(mean.pca.family.net2$family))]
                df[is.na(rowSums(df[, 2:3])), 2:3] <- 0
                if(sum(df[, 2:3] != 0) > 0){
                    res.df <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                                     correl = cor(df$net1, df$net2, method = "spearman")
                                     )
                    res <- rbind(res, res.df)
                } else {
                    res.df <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                                     correl = 0
                                     )
                    res <- rbind(res, res.df)
                }
            } else {
                net1.full$genus <- sapply(net1.full$species, function(x){strsplit(as.character(x), split = "_")[[1]][1]})
                net2.full$genus <- sapply(net2.full$species, function(x){strsplit(as.character(x), split = "_")[[1]][1]})
                mean.pca.genus.net1 <- setNames(aggregate(net1.full$pca, by = list(net1.full$genus), FUN = mean), c("genus", "pca"))
                mean.pca.genus.net2 <- setNames(aggregate(net2.full$pca, by = list(net2.full$genus), FUN = mean), c("genus", "pca"))
                df <- data.frame(genus = sort(unique(c(as.character(mean.pca.genus.net1$genus), as.character(mean.pca.genus.net2$genus)))),
                                 net1 = NA,
                                 net2 = NA)
                df$net1 <- mean.pca.genus.net1$pca[match(as.character(df$genus), as.character(mean.pca.genus.net1$genus))]
                df$net2 <- mean.pca.genus.net2$pca[match(as.character(df$genus), as.character(mean.pca.genus.net2$genus))]
                df[is.na(rowSums(df[, 2:3])), 2:3] <- 0
                if(sum(df[, 2:3] != 0) > 0){
                    res.df <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                                     correl = cor(df$net1, df$net2, method = "spearman")
                                     )
                    res <- rbind(res, res.df)
                } else {
                    res.df <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                                     correl = 0
                                     )
                    res <- rbind(res, res.df)
                }
            }
        }
    }
    res <- res[-1,]
    rand.dist.correl <-
        ggplot(data = res) +
        geom_point(aes(x = dist, y = correl), size = 2) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
        geom_smooth(aes(x = dist, y = correl), se = FALSE) +
        ylim(-1, 1)
    rand.dist.correl.build <- ggplot_build(rand.dist.correl)$data[[3]][, c("x", "y")]
    return(list(sims = res, data.plot = rand.dist.correl.build))
}


## Randomizing networks instead of species

rand.geog.nets.pernet.single0 <- function(x, level = "family"){
    res <- data.frame(dist = NA, correl = NA)
    metadata <- net.metadata
    metadata[, c("long_dec", "lat_dec")] <- metadata[sample(1:nrow(metadata)), c("long_dec", "lat_dec")]
    for(i in 1:28){
        for(j in (i + 1):29){
            net1.full <- fulldata.net.bin[fulldata.net.bin$id == unique(fulldata.net.bin$id)[i],]
            net2.full <- fulldata.net.bin[fulldata.net.bin$id == unique(fulldata.net.bin$id)[j],]
            if(level == "family"){
                mean.pca.family.net1 <- setNames(aggregate(net1.full$pca, by = list(net1.full$family), FUN = mean), c("family", "pca"))
                mean.pca.family.net2 <- setNames(aggregate(net2.full$pca, by = list(net2.full$family), FUN = mean), c("family", "pca"))
                df <- data.frame(family = sort(unique(c(as.character(mean.pca.family.net1$family), as.character(mean.pca.family.net2$family)))),
                                 net1 = NA,
                                 net2 = NA)
                df$net1 <- mean.pca.family.net1$pca[match(as.character(df$family), as.character(mean.pca.family.net1$family))]
                df$net2 <- mean.pca.family.net2$pca[match(as.character(df$family), as.character(mean.pca.family.net2$family))]
                df$net1[is.na(df$net1)] <- 0
                df$net2[is.na(df$net2)] <- 0
                res.df <- data.frame(dist = distm(x = c(metadata$long_dec[i], metadata$lat_dec[i]), y = c(metadata$long_dec[j], metadata$lat_dec[j]), fun = distGeo),
                                     correl = cor(df$net1, df$net2, method = "spearman")
                                     )
                res <- rbind(res, res.df)
            } else {
                net1.full$genus <- sapply(net1.full$species, function(x){strsplit(as.character(x), split = "_")[[1]][1]})
                net2.full$genus <- sapply(net2.full$species, function(x){strsplit(as.character(x), split = "_")[[1]][1]})
                mean.pca.genus.net1 <- setNames(aggregate(net1.full$pca, by = list(net1.full$genus), FUN = mean), c("genus", "pca"))
                mean.pca.genus.net2 <- setNames(aggregate(net2.full$pca, by = list(net2.full$genus), FUN = mean), c("genus", "pca"))
                df <- data.frame(genus = sort(unique(c(as.character(mean.pca.genus.net1$genus), as.character(mean.pca.genus.net2$genus)))),
                                 net1 = NA,
                                 net2 = NA)
                df$net1 <- mean.pca.genus.net1$pca[match(as.character(df$genus), as.character(mean.pca.genus.net1$genus))]
                df$net2 <- mean.pca.genus.net2$pca[match(as.character(df$genus), as.character(mean.pca.genus.net2$genus))]
                df$net1[is.na(df$net1)] <- 0
                df$net2[is.na(df$net2)] <- 0
                res.df <- data.frame(dist = distm(x = c(metadata$long_dec[i], metadata$lat_dec[i]), y = c(metadata$long_dec[j], metadata$lat_dec[j]), fun = distGeo),
                                     correl = cor(df$net1, df$net2, method = "spearman")
                                     )
                res <- rbind(res, res.df)
            }
        }
    }
    res <- res[-1,]
    rand.dist.correl <-
        ggplot(data = res) +
        geom_point(aes(x = dist, y = correl), size = 2) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
        geom_smooth(aes(x = dist, y = correl), se = FALSE) +
        ylim(-1, 1)
    rand.dist.correl.build.pernet <- ggplot_build(rand.dist.correl)$data[[3]][, c("x", "y")]
    return(list(sims = res, data.plot = rand.dist.correl.build.pernet))
}


rand.geog.nets.pernet.full0 <- function(x, level = "family"){
    res <- data.frame(dist = NA, correl = NA)
    metadata <- net.metadata
    metadata[, c("long_dec", "lat_dec")] <- metadata[sample(1:nrow(metadata)), c("long_dec", "lat_dec")]
    for(i in 1:28){
        for(j in (i + 1):29){
            net.index <- sample(seq(1, 29), 2)
            net1.full <- fulldata.net.bin[fulldata.net.bin$id == unique(fulldata.net.bin$id)[i],]
            net2.full <- fulldata.net.bin[fulldata.net.bin$id == unique(fulldata.net.bin$id)[j],]
            if(level == "family"){
                mean.pca.family.net1 <- setNames(aggregate(net1.full$pca, by = list(net1.full$family), FUN = mean), c("family", "pca"))
                mean.pca.family.net2 <- setNames(aggregate(net2.full$pca, by = list(net2.full$family), FUN = mean), c("family", "pca"))
                df <- data.frame(family = sort(unique(c(as.character(mean.pca.family.net1$family), as.character(mean.pca.family.net2$family)))),
                                 net1 = NA,
                                 net2 = NA)
                df$net1 <- mean.pca.family.net1$pca[match(as.character(df$family), as.character(mean.pca.family.net1$family))]
                df$net2 <- mean.pca.family.net2$pca[match(as.character(df$family), as.character(mean.pca.family.net2$family))]
                df[is.na(rowSums(df[, 2:3])), 2:3] <- 0
                if(sum(df[, 2:3] != 0) > 0){
                    res.df <- data.frame(dist = distm(x = c(metadata$long_dec[i], metadata$lat_dec[i]), y = c(metadata$long_dec[j], metadata$lat_dec[j]), fun = distGeo),
                                     correl = cor(df$net1, df$net2, method = "spearman")
                                     )
                    res <- rbind(res, res.df)
                } else {
                    res.df <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                                     correl = 0
                                     )
                    res <- rbind(res, res.df)
                }
            } else {
                net1.full$genus <- sapply(net1.full$species, function(x){strsplit(as.character(x), split = "_")[[1]][1]})
                net2.full$genus <- sapply(net2.full$species, function(x){strsplit(as.character(x), split = "_")[[1]][1]})
                mean.pca.genus.net1 <- setNames(aggregate(net1.full$pca, by = list(net1.full$genus), FUN = mean), c("genus", "pca"))
                mean.pca.genus.net2 <- setNames(aggregate(net2.full$pca, by = list(net2.full$genus), FUN = mean), c("genus", "pca"))
                df <- data.frame(genus = sort(unique(c(as.character(mean.pca.genus.net1$genus), as.character(mean.pca.genus.net2$genus)))),
                                 net1 = NA,
                                 net2 = NA)
                df$net1 <- mean.pca.genus.net1$pca[match(as.character(df$genus), as.character(mean.pca.genus.net1$genus))]
                df$net2 <- mean.pca.genus.net2$pca[match(as.character(df$genus), as.character(mean.pca.genus.net2$genus))]
                df[is.na(rowSums(df[, 2:3])), 2:3] <- 0
                if(sum(df[, 2:3] != 0) > 0){
                    res.df <- data.frame(dist = distm(x = c(metadata$long_dec[i], metadata$lat_dec[i]), y = c(metadata$long_dec[j], metadata$lat_dec[j]), fun = distGeo),
                                     correl = cor(df$net1, df$net2, method = "spearman")
                                     )
                    res <- rbind(res, res.df)
                } else {
                    res.df <- data.frame(dist = distm(x = c(net.metadata$long_dec[i], net.metadata$lat_dec[i]), y = c(net.metadata$long_dec[j], net.metadata$lat_dec[j]), fun = distGeo),
                                     correl = 0
                                     )
                    res <- rbind(res, res.df)
                    }
            }
        }
    }
    res <- res[-1,]
    rand.dist.correl <-
        ggplot(data = res) +
        geom_point(aes(x = dist, y = correl), size = 2) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
        geom_smooth(aes(x = dist, y = correl), se = FALSE) +
        ylim(-1, 1)
    rand.dist.correl.build.pernet <- ggplot_build(rand.dist.correl)$data[[3]][, c("x", "y")]
    return(list(sims = res, data.plot = rand.dist.correl.build.pernet))
}





registerDoMC(6)
correl.family.within.single0 <- foreach(k = 1:1000) %dopar% {
    print(k)
    rand.geog.nets.single0(fulldata.net.bin)
    }

correl.family.within.full0 <- foreach(k = 1:1000) %dopar% {
    print(k)
    rand.geog.nets.full0(fulldata.net.bin)
    }

correl.genus.within.single0 <- foreach(k = 1:1000) %dopar% {
    print(k)
    rand.geog.nets.single0(fulldata.net.bin, level = "genus")
    }

correl.genus.within.full0 <- foreach(k = 1:1000) %dopar% {
    print(k)
    rand.geog.nets.full0(fulldata.net.bin, level = "genus")
    }


correl.family.between.single0 <- foreach(k = 1:1000) %dopar% {
    print(k)
    rand.geog.nets.pernet.single0(fulldata.net.bin)
    }

correl.family.between.full0 <- foreach(k = 1:1000) %dopar% {
    print(k)
    rand.geog.nets.pernet.full0(fulldata.net.bin)
    }

correl.genus.between.single0 <- foreach(k = 1:1000) %dopar% {
    print(k)
    rand.geog.nets.pernet.single0(fulldata.net.bin, level = "genus")
    }

correl.genus.between.full0 <- foreach(k = 1:1000) %dopar% {
    print(k)
    rand.geog.nets.pernet.full0(fulldata.net.bin, level = "genus")
    }


null.smooth.family.within.single0 <- ldply(correl.family.within.single0, function(x){bind_rows(x[[2]])})
null.smooth.family.within.single0$sim <- rep(1:1000, each = 80)

null.smooth.family.within.full0 <- ldply(correl.family.within.full0, function(x){bind_rows(x[[2]])})
null.smooth.family.within.full0$sim <- rep(1:1000, each = 80)

null.smooth.genus.within.single0 <- ldply(correl.genus.within.single0, function(x){bind_rows(x[[2]])})
null.smooth.genus.within.single0$sim <- rep(1:1000, each = 80)

null.smooth.genus.within.full0 <- ldply(correl.genus.within.full0, function(x){bind_rows(x[[2]])})
null.smooth.genus.within.full0$sim <- rep(1:1000, each = 80)



null.smooth.family.between.single0 <- ldply(correl.family.between.single0, function(x){bind_rows(x[[2]])})
null.smooth.family.between.single0$sim <- rep(1:1000, each = 80)

null.smooth.family.between.full0 <- ldply(correl.family.between.full0, function(x){bind_rows(x[[2]])})
null.smooth.family.between.full0$sim <- rep(1:1000, each = 80)

null.smooth.genus.between.single0 <- ldply(correl.genus.between.single0, function(x){bind_rows(x[[2]])})
null.smooth.genus.between.single0$sim <- rep(1:1000, each = 80)

null.smooth.genus.between.full0 <- ldply(correl.genus.between.full0, function(x){bind_rows(x[[2]])})
null.smooth.genus.between.full0$sim <- rep(1:1000, each = 80)


save(null.smooth.family.within.single0,
     null.smooth.family.within.full0,
     null.smooth.genus.within.single0,
     null.smooth.genus.within.full0,
     null.smooth.family.between.single0,
     null.smooth.family.between.full0,
     null.smooth.genus.between.single0,
     null.smooth.genus.between.full0,
     file = "../output/null_correl_geog_dist_Mar2021.RData")
