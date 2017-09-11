####
#### R code for Ushio et al.
#### "Fluctuating interaction network and time-varying stability of a natural fish community"
#### No.4 Calculation of the dynamic stability
####

# Load config
source("config.R")

# Make output directory
kOutDir.03 <- "03_Network_out"
kOutDir.04 <- "04_Stability_out"
dir.create(kOutDir.04, showWarnings = FALSE)
# Save output of analysis in .RData file
kSaveOutput <- T

# Load previous work space
out.wd3 <-
  file.path(kOutDir.03,
            paste0("03_Network", config$kBestE.RangeStr, "_out.RData"))
load(out.wd3)
out.wd4 <-
  file.path(kOutDir.04,
            paste0("04_Stability", config$kBestE.RangeStr, "_out.RData"))

# Load packages
library(rEDM)
library(tseriesChaos)

# Load functions
source("functions/Functions_1_HelperFuncs.R")

# Extract species with no detected causal interactions and no links
new.nums <- new.nums.only.int
new.nums$xmap_from <- as.numeric(new.nums$xmap_from)
new.nums$xmap_to   <- as.numeric(new.nums$xmap_to)
check.link.spp   <- new.nums[new.nums$xmap_from != new.nums$xmap_to,]
linked.cause.spp <- unique(check.link.spp$xmap_from)
linked.all.spp   <-
  sort(unique(c(
    check.link.spp$xmap_from, check.link.spp$xmap_to
  )))
all.spp      <- 1:15
no.cause.spp <- all.spp[-linked.cause.spp]
no.link.spp  <- all.spp[-linked.all.spp]


# Make smap matrix (excluding intercept and lag)
smapc.no.const <- data.frame(tmp = rep(NA, nrow(smapc1)))
for (i in 1:ncol(smapc1)) {
  if (substr(colnames(smapc1)[i], 1, 9) == "const_for" |
      substr(colnames(smapc1)[i], 12, 14) == "lag") {

  } else{
    smapc.no.const <- cbind(smapc.no.const, smapc1[, i])
    colnames(smapc.no.const)[ncol(smapc.no.const)] <-
      colnames(smapc1)[i]
  }
}
smapc.noC.noLag <- smapc.no.const[,-1]

# Make smap coefficients with lags
smapc.lag <- data.frame(tmp = rep(NA, nrow(smapc1)))
for (i in 1:ncol(smapc1)) {
  if (substr(colnames(smapc1)[i], 12, 14) == "lag") {
    smapc.lag <- cbind(smapc.lag, smapc1[, i])
    colnames(smapc.lag)[ncol(smapc.lag)] <- colnames(smapc1)[i]
  }
}
smapc.Lag <- smapc.lag[,-1]


# Calculation of the dynamic stability
eigen.all <- data.frame(matrix(NA, ncol = 1))
colnames(eigen.all) <- "max_s"
for (i in 1:nrow(smapc.noC.noLag)) {
  j0 <- MakeSmapcMatrix2(new.nums, smapc.noC.noLag[i,])

  if (any(is.na(j0))) {
    eigen.all <- rbind(eigen.all, NA)
  } else{
    j0    <- j0[linked.all.spp, linked.all.spp]
    j.top <- j0
    for (j in 1:(max(config$kBestE.Range) - 1)) {
      j.temp <- ComputeJi(new.nums2, j, i)[linked.all.spp, linked.all.spp]
      j.top  <- cbind(j.top, j.temp)
    }

    cr.len  <- length(d.name[linked.all.spp])
    dim.is  <- max(config$kBestE.Range) - 1
    unities <- diag(cr.len * dim.is)
    zeros   <- matrix(0, nrow = cr.len * dim.is, ncol = cr.len)
    j.all   <- rbind(j.top, cbind(unities, zeros))

    # Calculating eigenvalues
    ev <- eigen(j.all)$values
    re.ev  <- Re(ev)
    Max_re <- max(abs(re.ev))
    max.ev <- abs(Re(ev[abs(Re(ev)) == Max_re]))

    if (length(max.ev) == 1) eigen_tmp <- max.ev
    if (length(max.ev) != 1) eigen_tmp <- max.ev[1]

    names(eigen_tmp) <- "max_s"
    eigen.all <- rbind(eigen.all, eigen_tmp)
  }
}
eigen.all <- eigen.all[-1,]

# Show the stability figure
plot(
  eigen.all ~ c(1:length(eigen.all)),
  type = "l",
  xlab = "Time",
  ylab = "Dynamic stability",
  axes = F,
  main = "Dynamic stability of Maizuru Fish Community"
)
box()
axis(2, las = 2); axis(1, las = 1)
abline(h = 1, lty = 2)
abline(v = (1:12) * 24,
       col = "gray50",
       lwd = 1)

# Save the result of the analysis
if (kSaveOutput) {
  # save configuration as config.04 (and not config)
  config.04 <- config
  rm(config)

  save.image(out.wd4)
}
