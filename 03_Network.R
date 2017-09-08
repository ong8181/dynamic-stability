####
#### R code for Ushio et al.
#### "Fluctuating interaction network and time-varying stability of a natural fish community"
#### No.3 Illustration of the interaction network
####

# Load config
source("config.R")

# Make output directory
kOutDir.02 <- "02_SmapCoef_out"
kOutDir.03 <- "03_Network_out"
dir.create(kOutDir.03, showWarnings = FALSE)

# Load previous work space
ws.out2 <-
  file.path(kOutDir.02,
            paste0("02_SmapCoef", config$kBestE.RangeStr, "_out.RData"))
load(ws.out2)
out.wd3 <-
  file.path(kOutDir.03,
            paste0("03_Network", config$kBestE.RangeStr, "_out.RData"))

# Save output of analysis in .RData file
kSaveOutput <- T

# Load packages
library(rEDM)
library(tseriesChaos)

# load functions
source("functions/Functions_1_HelperFuncs.R")

# Packages to write interaction network
library(GGally)
library(intergraph)
library(network)
library(igraph)
library(ggplot2)
library(gridExtra)

d.name.num <- data.frame(cbind(d.name, 1:length(d.name)))
d.name.num$combined  <- NA
for (i in 1:length(d.name)) {
  d.name.num$combined[i]   <-
    sprintf("%s_%s", d.name.num[i, 2], d.name.num[i, 1])
}

colnames(smapc.mat1) <- d.name.num$d.name
rownames(smapc.mat1) <- d.name.num$d.name

# Exclude self-regulation terms
diag(smapc.mat1) <- 0
net1 <- network(smapc.mat1)

is.no.self <-
  new.nums.only.int[new.nums.only.int$xmap_from != new.nums.only.int$xmap_to,]
n1 <- PlotNetwork(net1, is.no.self$tp1_m)
n1 <- n1 + labs(title = "Maizuru Fish community network")

# Show network figure
print(n1)

# Save the result of the analysis
if (kSaveOutput) {
  # save configuration as config.03 (and not config)
  config.03 <- config
  rm(config)
  
  save.image(out.wd3)
}
