####
#### R code for Ushio et al.
#### "Fluctuating interaction network and time-varying stability of a natural fish community"
#### No.0 CCM for Maizuru Fish community
####

# Load config
source("config.R")

# Load functions
source("functions/Functions_1_HelperFuncs.R")
source("functions/Functions_3_TwinSurrogate.R")
source("functions/Functions_4_Twin95CIpara.R")

# Load packages
library("rEDM"); packageVersion("rEDM")
library("tseriesChaos")
library("pforeach")

# Make output directory
kOutDir <- "00_MaizuruCCM_out_recalculate"
dir.create(kOutDir, showWarnings = FALSE)

# Load data
# Biweekly-observed dominant species data
biw.data <- read.csv('data/Maizuru_dominant_sp.csv')
d.name   <- names(biw.data)[4:18]


# Do CCM and check significance based on Phase-Lock Twin surrogate
m <- nrow(biw.data)
lib.size <- c(11, 15, 24, 40, 80, 120, 160, 200, 240, m)

# Best E is used, criteria = MAE, result saved in one figure
bestE.dir <- file.path(kOutDir,
                       paste0("rho_values_bestE_byMAE_",
                              config$kBestE.RangeStr))
dir.create(bestE.dir, showWarning = F)

# Do CCM for all species pairs
for (i in 1:length(d.name)) {
  for (j in 1:length(d.name)) {
    t.begin <- proc.time()

    xmap.from.to.name <- sprintf("Xmap from %s to %s", d.name[i], d.name[j])

    x <- as.numeric(scale(biw.data[, d.name[i]]))
    y <- as.numeric(scale(biw.data[, d.name[j]]))
    e.x <- BestErEDM(x, E = config$kBestE.Range)
    e.y <- BestErEDM(y, E = config$kBestE.Range)

    block <- cbind(x, y)

    ccm.raw <-
      ccm(
        block,
        E         = e.x,
        lib_sizes = lib.size,
        silent    = T
      )
    ccm.tmp <- ccm_means(ccm.raw)

    x.sur <- TwinSurrogate(x, e.x, 100, surrogate.option = "phase_lock")
    twin.sur <- ParSurCI(x, y, x.sur, E.range = config$kBestE.Range)

    twin.sur$rho$rho <- ccm.tmp$rho
    twin.sur$rho$u_rho <- ccm.tmp$rho - twin.sur$rho$upper95
    twin.sur$rho$d_rho <- ccm.tmp$rho[nrow(ccm.tmp)] - ccm.tmp$rho[1]
    ccm.tmp <- cbind(ccm.tmp, twin.sur$rho)

    ccm.list <- list(
      ccm_values  = ccm.tmp,
      max_rho     = max(ccm.tmp$rho),
      ter_rho     = rev(ccm.tmp$rho)[1],
      max_min_rho = (max(ccm.tmp$rho) - min(ccm.tmp$rho)),
      ini_max_rho = (max(ccm.tmp$rho) - ccm.tmp$rho[1]),
      d_rho       = (rev(ccm.tmp$rho)[1] - ccm.tmp$rho[1]),
      u_rho_term  = rev(twin.sur$rho$u_rho)[1],
      xmap_from = i,
      xmap_to   = j,
      xmap_from_to = xmap.from.to.name
    )

    ccm.res.file <- file.path(bestE.dir, sprintf("rhos_from_%g%s_to_%g%s.csv", i, d.name[i], j, d.name[j]))

    write.csv(ccm.list, ccm.res.file, row.names = F)
    time.used <- proc.time() - t.begin
    cat("CCM, Xmap from", i, "to", j, "completed:", time.used[3],"sec\n")
  }
}


# Save result
ws.out <-
  file.path(kOutDir,
            paste0(
              "01_MaizuruCCM_out_",
              config$kBestE.RangeStr,
              ".RData"
            ))

# save configuration as config.01 (and not config)
config.00 <- config
rm(config)

save.image(file = ws.out)
