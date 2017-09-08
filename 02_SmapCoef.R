####
#### R code for Ushio et al.
#### "Fluctuating interaction network and time-varying stability of a natural fish community"
#### No.2 Caluculate S-map Coefficients
####

# Load config
source("config.R")

# Save output of analysis in .RData file
kSaveOutput <- T

# Make output directory
kOutDir.01 <- "01_ExtractCausality_out"
kOutDir.02 <- "02_SmapCoef_out"
dir.create(kOutDir.02, showWarnings = FALSE)

# Load packages
library(rEDM)
library(tseriesChaos)
library(pforeach)

# Load functions
source("functions/Functions_1_HelperFuncs.R")
source("functions/Functions_2_SmapCoefPara.R")

# Load previous work space
ws.out1 <- file.path(kOutDir.01, paste0("01_ExtractCausality", config$kBestE.RangeStr, "_out.RData"))
load(ws.out1)
ws.out2 <- file.path(kOutDir.02, paste0("02_SmapCoef", config$kBestE.RangeStr, "_out.RData"))

# Rearreange the order of xmap_from and xmap_to
new.nums <- num.for.smapc
for (i in 1:length(d.name)) {
  num.tmp <- num.for.smapc[num.for.smapc$xmap_from == i,]
  tmp1 <- num.tmp[num.tmp$xmap_from[1] == num.tmp$xmap_to,]
  tmp2 <- num.tmp[num.tmp$xmap_from[1] != num.tmp$xmap_to,]
  new.nums[as.numeric(rownames(num.tmp)),] <- rbind(tmp1, tmp2)
}

# Quantify interaction strengths among the fish community
smapc.res <- SmapCFunc(new.nums, smapc.tp = 1, stats.output = T,
                       embedding = "best_E", original.data = biw.data)
smapc.tp1 <- smapc.res$coefs
smapc1 <- smapc.tp1[,-1]
smapc.tp1.m0 <- as.data.frame(matrix(colMeans(smapc1, na.rm = T), ncol = 1))
rownames(smapc.tp1.m0) <- colnames(smapc1)

smapc.tp1.m0$include <- smapc.tp1.m0[, 1]
for (i in 1:nrow(smapc.tp1.m0)) {
  if (substr(rownames(smapc.tp1.m0)[i], 1, 9) == "const_for") {
    smapc.tp1.m0$include[i] <- NA
  }
}
smapc.tp1.m <- na.omit(smapc.tp1.m0)


# Add lag effects
new.nums2 <- data.frame(NULL)
for (i in 1:length(d.name)) {
    realz.E <- nrow(new.nums[new.nums$xmap_from == i,])
    if (realz.E > 0) {
      embed.E <- unique(new.nums[new.nums$xmap_from == i, 'best_E'])
      add.E <- embed.E - realz.E
      if (add.E > 0) {
        xmap.from.name <- sprintf("lag%s", 1:add.E)
        new.nums.tmp1 <-
          data.frame(xmap_from = xmap.from.name,
                     xmap_to = i,
                     best_E = embed.E)
        new.nums.tmp2 <- new.nums[new.nums$xmap_from == i,]
        new.nums.tmp3 <- rbind(new.nums.tmp2, new.nums.tmp1)
      } else {
        new.nums.tmp3 <- new.nums[new.nums$xmap_from == i,]
      }
      new.nums2 <- rbind(new.nums2, new.nums.tmp3)
    }
}

rownames(new.nums2) <- rownames(smapc.tp1.m)
new.nums2$tp1_m <- smapc.tp1.m[, 1]
new.nums2$tp1_m_name <- rownames(smapc.tp1.m)

is.lag <- substr(new.nums2$tp1_m_name, 12, 14) == "lag"
new.nums.only.int <- new.nums2[!is.lag,]
smapc.mat1 <- MakeSmapcMatrix(new.nums.only.int)

# Save the result of the analysis
if (kSaveOutput) {
  # save configuration as config.02 (and not config)
  config.02 <- config
  rm(config)

  save.image(ws.out2)
}
