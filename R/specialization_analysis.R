library("tidyverse")
library("taxize")
library("picante")
source("./R/calculating_centrality.R")

plant.families <- data.frame(
    sp = as.character(setNames(unlist(sapply(nets, names)), NULL)),
    stringsAsFactors = FALSE
)

plant.families$genus <- sapply(strsplit(plant.families$sp, split = ".", fixed = TRUE), function(x){x[[1]][1]})

for(i in 1:nrow(plant.families)){
    print(i)
    plant.families$family[i] <- ipni_search(genus = plant.families$genus[i])$family[1]
}

plant.families$family[is.na(plant.families$family)] <-
    c(NA, "Primulaceae", NA,
      "Primulaceae", "Bromeliaceae", NA,
      "Marcgraviaceae", NA, NA,
      NA, "Loranthaceae", "Loranthaceae",
      NA, "Lauraceae", NA,
      NA, NA, "Phyllanthaceae",
      "Moraceae", NA, NA,
      "Melastomataceae", NA, NA,
      NA, NA, NA,
      NA, NA, NA,
      NA, NA, NA,
      NA, NA, "Anacardiaceae",
      "Magnoliaceae", "Melastomataceae", "Araceae",
      "Rhamnaceae", "Melastomataceae", "Celastraceae",
      "Rubiaceae", NA, "Cucurbitaceae",
      "Melastomataceae", "Melastomataceae", "Melastomataceae",
      "Melastomataceae", "Melastomataceae", "Loranthaceae",
      "Loranthaceae", "Vitaceae", NA,
      NA, NA, NA,
      NA, NA, NA,
      NA, NA, NA,
      NA, NA, NA,
      NA, NA, NA,
      NA, NA, NA,
      NA, NA, NA,
      NA, NA, NA,
      NA, "Rosaceae", NA,
      NA, NA, NA,
      NA, NA, "Bromeliaceae",
      "Bromeliaceae", "Rubiaceae", "Loranthaceae",
      NA, "Lauraceae", "Poaceae",
      "Smilacaceae")

plant.tree <- read.tree("./data/plant_family_tree.tre")
plant.tree$node.label <- NULL

for(i in 1:nrow(fulldata.net.bin)){
    print(i)
    net.id <- match(fulldata.net.bin$id[i], unique(fulldata.net.bin$id))
    inter.sp <- nets[[net.id]][match(fulldata.net.bin$species[i], rownames(nets[[net.id]])),]
    families <- unique(plant.families$family[match(names(inter.sp)[which(inter.sp > 0)], plant.families$sp)])
    if(!all(is.na(families))){
        fulldata.net.bin$n.plant.families[i] <- length(na.omit(families))
        tr.temp <- drop.tip(plant.tree, tip = plant.tree$tip.label[is.na(match(plant.tree$tip.label, na.omit(families)))])
        if(!is.null(tr.temp)){
            fulldata.net.bin$phylodiv[i] <- ifelse(Ntip(tr.temp) == 1, 0, sum(tr.temp$edge.length))
        } else {
            fulldata.net.bin$phylodiv[i] <- NA
        }
    } else {
        fulldata.net.bin$phylodiv[i] <- NA
    }
}


write.table(fulldata.net.bin, file = "./data/interaction_per_taxonomy.csv", quote = FALSE, row.names = FALSE, sep = ",")

