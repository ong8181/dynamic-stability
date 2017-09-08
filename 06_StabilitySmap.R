####
#### R code for Ushio et al.
#### "Fluctuating interaction network and time-varying stability of a natural fish community"
#### No.6 Calculate influences from network properties to the dynamic stability
####

# Load config
source("config.R")

# Make output directory
kOutDir.05 <- "05_StabilityCCM_out"
kOutDir.06 <- "06_StabilitySmap_out"
dir.create(kOutDir.06, showWarnings = FALSE)

# Load previous work space
out.wd5 <-
  file.path(kOutDir.05,
            paste0("05_StabilityCCM_out_", config$kBestE.RangeStr, ".RData"))
load(out.wd5)
out.wd6 <-
  file.path(
    kOutDir.06,
    paste0(
      "06_StabilitySmap_out_",
      config$kBestE.RangeStr,
      ".RData"
    )
  )


# Load packages
library("rEDM")
library("tseriesChaos")

# Load functin
source("functions/Functions_1_HelperFuncs.R")

# Perform multivariate S-map
smapc.stability <- StabilitySmap(stability, cbind(int_mean_abs, weak_index, simpson))

# Plot results
int.mean.smapc <- smapc.stability[[1]]$smap_coefficients[,2]
weak.smapc <- smapc.stability[[1]]$smap_coefficients[,3]
simpson.smapc <- smapc.stability[[1]]$smap_coefficients[,4]

boxplot(smapc.stability[[1]]$smap_coefficients[,2:4],
        names = c("Mean IS", "Weak index", "Simpson's index"),
        ylab = "Influence to the stability",
        main = "Influences of network properties on the dynamic stability")
abline(h = 0, lty = 2)


# Save the result of the analysis
if (kSaveOutput) {
  # save configuration as config.05 (and not config)
  config.06 <- config
  rm(config)
  
  save.image(out.wd6)
}
