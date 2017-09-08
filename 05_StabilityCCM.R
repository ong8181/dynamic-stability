####
#### R code for Ushio et al.
#### "Fluctuating interaction network and time-varying stability of a natural fish community"
#### No.5 CCM between the dynamic stability and other variables
####

# Load config
source("config.R")

# Make output directory
kOutDir.04 <- "04_Stability_out"
kOutDir.05 <- "05_StabilityCCM_out"
dir.create(kOutDir.05, showWarnings = FALSE)

# Load previous work space
out.wd4 <-
  file.path(kOutDir.04,
            paste0("04_Stability", config$kBestE.RangeStr, "_out.RData"))
load(out.wd4)
out.wd5 <-
  file.path(
    kOutDir.05,
    paste0(
      "05_StabilityCCM_out_",
      config$kBestE.RangeStr,
      ".RData"
    )
  )


# Load packages
library("rEDM")
library("tseriesChaos")
library("pforeach")

# Calculations of community indices
self.smap <-
  sprintf('smapc_from_%s_to_%s', 1:length(d.name), 1:length(d.name))
no.self <- is.na(match(colnames(smapc.noC.noLag), self.smap))
smapc.noC.noLag.noSelf <- smapc.noC.noLag[, no.self]
interaction.mat <- smapc.noC.noLag.noSelf
ds <- cbind(biw.data, interaction.mat, eigen_all = eigen.all)

# From biw.data
sp.col <- 4:(ncol(biw.data) - 3)
sp.end <- ncol(biw.data)
int.names <-
  new.nums[new.nums$xmap_from != new.nums$xmap_to, 'tp1_m_name']
int.end   <- sp.end + ncol(interaction.mat)
eigen.col <- int.end + 1
d.dominant  <- ds[, d.name[linked.all.spp]]


# Calculate community indices
ds$total_dom <- rowSums(d.dominant)
prop.tmp <- d.dominant / ds$total_dom
ds$simp_dom <- apply(prop.tmp, 1, function(x) 1 - sum(x^2))
int.cols <- which(substr(colnames(ds), 1, 11) == "smapc_from_")

ds <- na.omit(ds)
ds$int_mean_abs <- apply(abs(ds[,int.cols]), 1, function(x) mean(x, na.rm=T))
ds$int_median_abs <- apply(abs(ds[,int.cols]), 1, function(x) median(x, na.rm=T))
ds$max_int_abs <- apply(abs(ds[,int.cols]), 1, function(x) max(x, na.rm=T))
ds$weak_index <- ds$int_median_abs/ds$max_int_abs

op <- par(mfrow = c(3,1))
plot(ds$int_mean_abs, type = "l", xlab = "Time index", ylab = "Mean IS")
plot(ds$weak_index, type = "l", xlab = "Time index", ylab = "Med/Max IS")
plot(ds$simp_dom, type = "l", xlab = "Time index", ylab = "Simpson's index")
par(op)


#### CCM between dynamic stability and interaction network properties
# Load functions
source("functions/Functions_1_HelperFuncs.R")
source("functions/Functions_3_TwinSurrogate.R")
source("functions/Functions_4_Twin95CIpara.R")

# The number of surrogate
kNN.S <- 100
# Save output of analysis in .RData file
kSaveOutput <- T

# Rename object
stability <- ds$eigen_all
weak_index <- ds$weak_index
int_mean_abs <- ds$int_mean_abs
simpson <- ds$simp_dom

int_res <- CcmStability(stability, int_mean_abs)
weak_res <- CcmStability(stability, weak_index)
simp_res <- CcmStability(stability, simpson)

op <- par(mfrow = c(1,3))
PlotStabilityCCM(int_res)
PlotStabilityCCM(weak_res)
PlotStabilityCCM(simp_res)
par(op)


# Save the result of the analysis
if (kSaveOutput) {
  # save configuration as config.05 (and not config)
  config.05 <- config
  rm(config)
  
  save.image(out.wd5)
}
