####
#### R code for Ushio et al.
#### "Fluctuating interaction network and time-varying stability of a natural fish community"
#### No.S1 Helper functions
####

# Check best embedding dimension
BestErEDM <- function(time.series,
                      E        = 1:30,
                      save.raw.data = F) {
  simplex.res <- simplex(time.series, E = E, exclusion_radius = 0)

  if (save.raw.data) {
    return(simplex.res)
  } else{
    bestE <- simplex.res[simplex.res$mae == min(simplex.res$mae), 'E']
    return(bestE)
  }
}

# For organizing output
CombineNames <- function(data) {
  data$from_to_names <- NA
  for (i in 1:nrow(data)) {
    data$from_to_names[i] <-
      sprintf("%s, from_%g_to_%g",
              data$xmap_from_to[i],
              data$from_num[i],
              data$to_num[i])
  }
  return(data)
}


# Network plot
PlotNetwork <- function(data, smapc.data, threshold = 0) {
  seg.sizes <- abs(smapc.data)
  smapc.col <- smapc.col.tmp <- smapc.data
  smapc.col[abs(smapc.col.tmp) < threshold]  <- "grey"
  smapc.col[smapc.col.tmp >= threshold] <- "royalblue"
  smapc.col[smapc.col.tmp <= -threshold]  <- "red3"

  return(
    ggnet(
      data,
      arrow.size    = 10,
      mode          = "circle",
      size          = 5,
      node.color    = "grey30",
      node.alpha    = 0.3,
      color         = "black",
      label.node    = T,
      label.size    = 4,
      segment.color = smapc.col,
      segment.alpha = 0.7,
      segment.size  = (seg.sizes * 3 + 1)
    )
  )
}


# matrix functions
MakeSmapcMatrix <- function(num_data) {
  mat.data <- matrix(rep(NA, length(d.name) ^ 2), ncol = length(d.name))
  for (i in 1:NROW(num_data)) {
    row.xmap.from <- as.numeric(num_data$xmap_from[i])
    col.xmap.to   <- as.numeric(num_data$xmap_to[i])
    # change row and col for network figure
    mat.data[col.xmap.to, row.xmap.from] <- num_data[i, 'tp1_m']
  }
  mat.data[is.na(mat.data)] <- 0
  colnames(mat.data) <- d.name
  rownames(mat.data) <- d.name
  return(mat.data)
}


MakeSmapcMatrix2 <- function(num.data.smapc, is.values) {
  if (any(is.na(is.values))) {
    return(NA)
  } else{
    nl   <- length(d.name)
    data <- matrix(rep(NA, nl ^ 2), ncol = nl)
    is.data <-
      cbind(num.data.smapc, as.numeric(is.values), names(is.values))
    colnames(is.data)[6] <- "temporal_IS"

    for (i in 1:NROW(is.data)) {
      row.xmap.from <- is.data$xmap_from[i]
      col.xmap.to   <- is.data$xmap_to[i]
      # change row and col for network figure
      data[row.xmap.from, col.xmap.to] <- is.data[i, 'temporal_IS']
    }

    data[is.na(data)] <- 0
    colnames(data) <- d.name
    rownames(data) <- d.name
    return(data)
  }
}


MakeSmapcMatrix3 <- function(num_data) {
  sub.name <- d.name[selected.spp]
  mat.data <- matrix(rep(NA, length(sub.name) ^ 2), ncol = length(sub.name))
  for (i in 1:NROW(num_data)) {
    row.xmap.from <- as.numeric(num_data$xmap_from[i])
    col.xmap.to   <- as.numeric(num_data$xmap_to[i])
    # Change row and col for network figure
    num.row <- match(row.xmap.from, sort(unique(as.numeric(num_data$xmap_from))))
    num.col <- match(col.xmap.to, sort(unique(as.numeric(num_data$xmap_to))))
    mat.data[num.col, num.row] <- num_data[i, 'tp1_m']
  }
  mat.data[is.na(mat.data)] <- 0
  colnames(mat.data) <- sub.name
  rownames(mat.data) <- sub.name
  return(mat.data)
}


ComputeJi <-
  function(smapc.data, lag.i, k, smapc.Lag.data = smapc.Lag) {
    ji <- matrix(rep(NA, length(d.name) ^ 2), ncol = length(d.name))
    for (i in 1:length(d.name)) {
      lag.name   <- sprintf("lag%s", lag.i)
      smapc.name <- smapc.data[smapc.data[, 'xmap_to'] == i &
                                 smapc.data[, 'xmap_from'] == lag.name,
                               "tp1_m_name"]
      smapc.Lag.tmp <- smapc.Lag.data[, smapc.name]
      if (length(smapc.Lag.tmp) < 1) {
        ji[i, i] <- 0
      }
      if (length(smapc.Lag.tmp) > 0) {
        ji[i, i] <- smapc.Lag.tmp[k]
      }
    }
    ji[is.na(ji)] <- 0
    return(ji)
  }

# Functions for CCM of the dynamic stability
CcmStability <- function(x, y){
  block0 <- cbind(x, y)
  block <- na.omit(block0)[,1:2]
  block <- apply(block, 2, function(x) as.numeric(scale(x)))

  e.x <- BestErEDM(block[,1], E = config$kBestE.Range)
  e.y <- BestErEDM(block[,2], E = config$kBestE.Range)

  m <- nrow(block)
  lib.size <- c(11, 15, 24, 40, 80, 120, 160, 200, 240, m)

  ccm.raw <-
    ccm(
      block,
      E         = e.x,
      lib_sizes = lib.size,
      silent    = T
    )

  ccm.tmp <- ccm_means(ccm.raw)
  x.sur <- TwinSurrogate(block[,1], e.x, kNN.S, surrogate.option = "phase_lock")
  twin.sur <- ParSurCI(block[,1], block[,2], x.sur, lib.parms = lib.size, E.range = config$kBestE.Range)
  result <- list(ccm_res = ccm.tmp, sur_res = twin.sur)

  return(result)
}

PlotStabilityCCM <- function(stability.ccm.res){
  plot(stability.ccm.res$ccm_res$lib_size,
       stability.ccm.res$ccm_res$rho,
       type = "l",
       col = "red3",
       xlab = "Library size",
       ylab = "CCM skill",
       ylim = c(0,1))
  polygon(c(stability.ccm.res$sur_res$rho$L,
            rev(stability.ccm.res$sur_res$rho$L)),
          c(stability.ccm.res$sur_res$rho$lower95,
            rev(stability.ccm.res$sur_res$rho$upper95)),
          border = F,
          col = rgb(0,0,0,0.2))
}


# Function for multivariate S-map for the dynamic stability
StabilitySmap <- function(effect.ts, cause.ts){
  y0 <- as.numeric(scale(effect.ts))
  if(length(cause.ts)>1) x0 <- apply(cause.ts, 2, function(x) as.numeric(scale(x)))
  if(length(cause.ts)<2) x0 <- as.numeric(scale(ds[,cause.ts]))
  block_0 <- cbind(y0, x0)
  Ey <- BestErEDM(y0, E = config$kBestE.Range)
  
  #### Add sufficient dimensions
  eE <- Ey - ncol(block_0) + 1
  NaN.mat <- matrix(rep(NaN, (eE - 1)*eE), ncol = eE)
  embed.ts <- rbind(NaN.mat, embed(y0, eE))
  block <- cbind(block_0, embed.ts[,2:eE])
  block <- na.omit(block)
  
  #### Do S-map
  # determine the best theta
  th.test <- block_lnlp(block, method="s-map", tp=1, theta=seq(0,10,by=0.1), silent=T, num_neighbors=0)
  best.th <- th.test[which.min(th.test$mae),'theta']
  
  #### Perform multivariate S-map to quantify interaction strength
  smapc.output <- block_lnlp(block,
                           method="s-map",
                           tp=1,
                           theta=best.th,
                           num_neighbors=0,
                           silent=T,
                           save_smap_coefficients=T)
  return(smapc.output)
}

