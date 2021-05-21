#### Calculating centrality for each species
library("plyr")
library("knitr")
library("rmarkdown")
library("igraph")
#library("qgraph")
library("tidyverse")
library("cowplot")
library("centiserve")
library("ape")
library("BAMMtools")
library("phangorn")
library("phytools")
library("raster")
library("sp")
library("rgdal")

## Importing networks and checking how many species of both birds and plants are in each.

nets <- lapply(as.list(paste0("../data/nets/", list.files("../data/nets/"))), read.csv, row.names = 1)
plyr::ldply(nets, dim)
names(nets) <- gsub(".csv", "", list.files("../data/nets/"))

## Excluding "problematic" networks (see documents for explanation)
nets <- nets[-c(2, 7, 10, 18, 27)]  # removed "BEEH", "CROM", "GEN_A", "LAMB", "SAPF"

## Removing species that are in separate modules in 6 networks that are not fully connected

## Removing from CACI - Tiaris_bicolor and Piper sp.
nets$CACI <- nets$CACI[-which(rownames(nets$CACI) == "Tiaris_bicolor"), -which(names(nets$CACI) == "Piper.sp.")]
## Removing from CACO - Euphonia_musica, Phoradendron.spp., Anthurium.scandens, Dendropemon.bicolor
nets$CACO <- nets$CACO[-which(rownames(nets$CACO) == "Euphonia_musica"), -match(c("Phoradendron.spp.", "Anthurium.scandens", "Dendropemon.bicolor"), names(nets$CACO))]
## Removing from CAFR - Vireo_latimeri, Casearia.sylvestris
nets$CAFR <- nets$CAFR[-which(rownames(nets$CAFR) == "Vireo_latimeri"), -which(names(nets$CAFR) == "Casearia.sylvestris")]
## Removing from MACK - Cracticus_cassicus + Heritiera.sp., Reinwardtoena_reinwardtsi + Osmoxylum.sp., Glycichaera_fallax + Poikilospermum.sp., Cicinnurus_magnificus + Philemon_buceroides + Sloanea.sp + Schefflera.spp. + Psychotria.sp., Melilestes_megarhynchus + Ficus.dammaropsis, Macropygia_phasianella + Symplocos.sp.
nets$MACK <- nets$MACK[-match(c("Cracticus_cassicus", "Reinwardtoena_reinwardtsi", "Glycichaera_fallax", "Cicinnurus_magnificus", "Philemon_buceroides", "Melilestes_megarhynchus", "Macropygia_phasianella"), rownames(nets$MACK)), -match(c("Heritiera.sp.", "Osmoxylum.sp.", "Poikilospermum.sp.", "Sloanea.sp.", "Schefflera.spp.", "Psychotriasp.", "Ficus.dammaropsis", "Symplocos.sp."), names(nets$MACK))]
## Removing from SAAV2 - Psarocolius_atrovirens + Dictyocaryum.lamarckianum
nets$SAAV_B <- nets$SAAV_B[-which(rownames(nets$SAAV_B) == "Psarocolius_atrovirens"), -which(names(nets$SAAV_B) == "Dictyocaryum.lamarckianum")]
## Removing from SARM - Dysithamnus_mentalis + SP23, Tangara_cayana + SP25, Euphonia_violacea + SP22
nets$SARM <- nets$SARM[-match(c("Dysithamnus_mentalis", "Tangara_cayana", "Euphonia_violacea"), rownames(nets$SARM)), -match(paste0("SP", c(22, 23, 25)), names(nets$SARM))]


### Network metrics

## Binarizing quantitative networks
bin.nets <- lapply(nets, function(x){x[x != 0] <- 1; return(x)})

## Calculating degree, closeness, Katz and page rank centrality

nets.graph.bin <- lapply(bin.nets, graph_from_incidence_matrix)
nets.pal <- RColorBrewer::brewer.pal(3, "Set2")[1:2]

## function to correctly align PCs

## foo <- function(x){
##     if(any(x$rotation[,1] < 0)){return(-x$x[,1])} else {return(x$x[,1])}
## }

df.pca.bin <- list()

for(i in 1:length(nets.graph.bin)){
    print(i)
    V(nets.graph.bin[[i]])$color <- ifelse(V(nets.graph.bin[[i]])$type == FALSE, nets.pal[2], nets.pal[1])
    V(nets.graph.bin[[i]])$degree <- igraph::degree(nets.graph.bin[[i]], mode = "all")
    V(nets.graph.bin[[i]])$closeness <- igraph::closeness(nets.graph.bin[[i]], mode = "all")
    V(nets.graph.bin[[i]])$katz <- centiserve::katzcent(nets.graph.bin[[i]], alpha = 0.05)
    df.pca.bin[[i]] <- data.frame(degree = V(nets.graph.bin[[i]])$degree, closeness = V(nets.graph.bin[[i]])$closeness, katz = V(nets.graph.bin[[i]])$katz)
    #V(nets.graph.bin[[i]])$pca <- foo(prcomp(df.pca.bin[[i]][, c("degree", "closeness", "katz")], scale = TRUE))
}


####################################################
## Creating table with all network data
####################################################

full.table.bin <- ldply(nets.graph.bin, function(x){data.frame(species = V(x)$name,
                                         degree = V(x)$degree,
                                         closeness = V(x)$closeness,
                                         katz = V(x)$katz,
                                         #pca = V(x)$pca,
                                         degree.z = scale(V(x)$degree),
                                         closeness.z = scale(V(x)$closeness),
                                         katz.z = scale(V(x)$katz)
                                         )})

full.table.bin <- full.table.bin[grep("_", full.table.bin$species),]
full.table.bin <- full.table.bin[-grep("\\.", full.table.bin$species),]

full.table.bin$pca <- prcomp(full.table.bin[, c("degree.z", "closeness.z", "katz.z")], scale = FALSE)$x[,1]

net.metadata <- read.csv("../data/network_data.csv")

full.table.bin$lat.clim <- net.metadata$Latitude[match(full.table.bin$.id, net.metadata$Network)]

### Importing lat/long data
net.metadata$lat_dec[net.metadata$lat_hem == "S"] <- -net.metadata$lat_dec[net.metadata$lat_hem == "S"]
net.metadata$long_dec[net.metadata$long_hem == "W"] <- -net.metadata$long_dec[net.metadata$long_hem == "W"]


### Obtaining climate data for each network

r <- raster::getData("worldclim", var="bio", res=10)

r <- r[[c(1, 4, 5, 6, 12, 15)]]
names(r) <- c("ann.mean.temp", "temp.seas", "max.warm.temp", "min.col.temp", "ann.prec", "prec.seas")

coords <- net.metadata[, c("long_dec", "lat_dec")]
names(coords) <- c("x", "y")

points <- SpatialPoints(coords)

values <- extract(r,points)

df <- cbind.data.frame(coordinates(points),values)
df[, 3:8] <- df[, 3:8]/10

xy <- df[, c(1,2)]

net.metadata <- cbind(net.metadata, df[, c(3, 4, 7, 8)])

### Including climate data on centrality table
fulldata.net.bin <- cbind(full.table.bin, net.metadata[match(full.table.bin$.id, as.character(net.metadata$Network)), -c(1, 6)])

names(fulldata.net.bin)[1] <- "id"

## write.table(fulldata.net.bin, file = "../data/full_centrality_bin.csv", sep = ",", quote = FALSE, row.names = FALSE)
