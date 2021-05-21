## Network data was organized, processed and exported using a separate script

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


## Removing species without rates

fulldata.net.bin <- fulldata.net.bin[!is.na(fulldata.net.bin$epsilon.eric0102),]

clim.pca <- prcomp(net.metadata[, 14:17], scale = TRUE)

fulldata.net.bin <- cbind(fulldata.net.bin, clim.pca$x[match(fulldata.net.bin$id, net.metadata$Network),])
names(fulldata.net.bin)[(ncol(fulldata.net.bin) - 3):ncol(fulldata.net.bin)] <- paste0("clim.pc", 1:4)

fulldata.net.bin <- cbind(fulldata.net.bin, setNames(lambda.split[match(fulldata.net.bin$species, lambda.split$species), -1], paste0("lambda.", names(lambda.split)[-1])))
fulldata.net.bin <- cbind(fulldata.net.bin, setNames(mu.split[match(fulldata.net.bin$species, mu.split$species), -1], paste0("mu.", names(mu.split)[-1])))
