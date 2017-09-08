####
#### R code for Ushio et al.
#### "Fluctuating interaction network and time-varying stability of a natural fish community"
#### No.1 Extract causality of Maizuru fish community
####

# Load config
source("config.R")

# Load data and CCM results
biw.data <- read.csv("data/Maizuru_dominant_sp.csv")
d.name   <- names(biw.data)[4:18]
d.bestE.0 <- read.csv("00_MaizuruCCM_out/Maizuru_ccm_res.csv")

# Make output directory
kOutDir.01 <- "01_ExtractCausality_out"
dir.create(kOutDir.01, showWarnings = FALSE)
ws.out1 <-
  file.path(kOutDir.01,
            paste0("01_ExtractCausality", config$kBestE.RangeStr, "_out.RData"))

# Save output of analysis in .RData file
kSaveOutput <- T
# Causality criteria
# Difference between upper limit of 95%CI and terminal rho
kURhoThreshold <- 0
# Threshold of delta rho (= terminal rho minus initial rho)
kDRhoThreshold <- 0.1

# Load packages
library(rEDM)
library(tseriesChaos)

# Load functions
source("functions/Functions_1_HelperFuncs.R")

# Load merged rho-value file
colnames(d.bestE.0)[14] <- 'from_num'
colnames(d.bestE.0)[15] <- 'to_num'

# setup data structure for human-readable output
d.bestE.1 <- CombineNames(d.bestE.0)

# Detect significant causal interactions
terminal.L <- max(d.bestE.1$ccm_values.L)
d.bestE    <- subset(d.bestE.1, ccm_values.L == terminal.L)

d.bestE.cause <- subset(d.bestE, ccm_values.u_rho > kURhoThreshold)
d.bestE.cause <- subset(d.bestE.cause, d_rho > kDRhoThreshold)
cause.pairs   <- unique(d.bestE.cause$from_to_names)

# Preparation for S-map coefficients analysis
num.for.smapc <- data.frame()
for (i in 1:length(cause.pairs)) {
  num.tmp <-
    d.bestE.cause[d.bestE.cause$from_to_names ==
                    cause.pairs[i], ][, c('from_num', 'to_num', 'ccm_values.E')]
  num.for.smapc <- rbind(num.for.smapc, num.tmp)
}

rownames(num.for.smapc) <- 1:nrow(num.for.smapc)
colnames(num.for.smapc) <- c('xmap_from', 'xmap_to', 'best_E')

# Save the result of the analysis
if (kSaveOutput) {
  # save configuration as config.01 (and not config)
  config.01 <- config
  rm(config)

  save.image(ws.out1)
}
