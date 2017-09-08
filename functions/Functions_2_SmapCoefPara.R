####
#### R code for Ushio et al.
#### "Fluctuating interaction network and time-varying stability of a natural fish community"
#### No.S2 Functions for calculating S-map coefficients
####

# Function that calculates s-map coefficient

# create custom lagged embedding
Embed2 <- function(ts, dim) {
  embed.c1 <- data.frame(lag0 = ts)
  if (dim > 1) {
    for (i in 1:(dim - 1)) {
      embed.c.tmp <- c(rep(NaN, i), ts[1:(length(ts) - i)])
      embed.c.tmp <- data.frame(matrix(embed.c.tmp, ncol = 1))
      embed.c1    <- cbind(embed.c1, embed.c.tmp)
      c.name      <- sprintf("lag%s", i)
      colnames(embed.c1)[ncol(embed.c1)] <- c.name
    }
  }
  return(embed.c1)
}

# SmapCoef for the converged pairs
SmapCFunc <- function(smapc.num.data,
                         smapc.tp     = 1,
                         embedding    = "min_multivariate",  # or "best_E"
                         stats.output = F,
                         original.data = biw.data) {
  # smapc.num.data <- new.nums
  sc.sum <- data.frame(tmp = rep(NA, NROW(original.data)))
  if (stats.output) {
    sc.stat <- data.frame(
      num_pred = rep(NA, NROW(smapc.num.data)),
      rho      = rep(NA, NROW(smapc.num.data)),
      mae      = rep(NA, NROW(smapc.num.data)),
      rmse     = rep(NA, NROW(smapc.num.data))
    )
  }

  smap.spp.examined <- 0

  for (i in 1:NROW(smapc.num.data)) {
    t.begin <- proc.time()
    # Preparation for fully multivariate S-map
    # Collect causal species
    effect.n   <- smapc.num.data[i, 'xmap_from']
    cause.main <- smapc.num.data[i, 'xmap_to']
    causal.spp <- smapc.num.data[smapc.num.data[, 'xmap_from'] == effect.n, 'xmap_to']
    smap.spp <- unique(c(causal.spp))

    if (smap.spp.examined != effect.n) {
      # Make a matrix for smap
      if (embedding == "best_E") {
        # add dimension
        eE <- smapc.num.data[i, 'best_E'] - length(smap.spp) + 1
        embed.ts <- Embed2(as.numeric(scale(original.data[, d.name[effect.n]])), eE)

        time.series.0 <- cbind(original.data[, d.name[smap.spp]])
        block.0 <-
          apply(time.series.0, 2, function(x)
            as.numeric(scale(x)))
        if (eE > 1) {
          block <- cbind(block.0, embed.ts[, 2:eE])
        }
        if (eE < 2) {
          block <- block.0
        }
      } else if (embedding == "min_multivariate") {
        time.series.0 <- cbind(original.data[, d.name[smap.spp]])
        block <-
          apply(time.series.0, 2, function(x)
            as.numeric(scale(x)))
      }

      # Determine the best theta
      theta.examined <- seq(0, 10, by = 0.1)
      th.test <-
        pforeach(
          i      = theta.examined,
          .c     = rbind,
          .cores = config$kMaxCore,
          .seed  = config$kRndSeed
        )({
          th.test0 <- block_lnlp(block, method = "s-map", tp = smapc.tp,
                                theta  = i, silent = T, num_neighbors = 0)
        })

      best.th <- th.test[th.test$mae == min(th.test$mae), 'theta']

      # Perform multivariate S-map to quantify interaction strengths
      smapc.res <- block_lnlp(block, method = "s-map", tp = smapc.tp,
                              theta  = best.th, num_neighbors = 0, silent = T,
                              save_smap_coefficients = T)
      smapc.tmp <- smapc.res[[1]]$smap_coefficients
      smapc.tmp <- as.data.frame(smapc.tmp)

      # Column names
      if (embedding == "best_E") {
        if (eE > 1) {
          smapc.col.names <- sprintf("smapc_from_%s_to_%s",
                                     c(smap.spp,
                                       colnames(embed.ts)[2:ncol(embed.ts)]),
                                     effect.n)
        } else if (eE < 2) {
          smapc.col.names <-
            sprintf("smapc_from_%s_to_%s", smap.spp, effect.n)
        }
      } else if (embedding == "min_multivariate") {
        smapc.col.names <-
          sprintf("smapc_from_%s_to_%s", smap.spp, effect.n)
      }
      smapc.col.names <-
        c(smapc.col.names, sprintf("const_for_%s", effect.n))
      colnames(smapc.tmp) <-  smapc.col.names
      sc.sum <- cbind(sc.sum, smapc.tmp)

      smap.spp.examined <- effect.n

      if (stats.output) {
        sc.stat[i,] <- smapc.res[[1]]$stats[, 1:4]
      }
    }

    t.used <- proc.time() - t.begin
    cat("Process ", i, "/", NROW(smapc.num.data), "completed:", t.used[3], "sec\n")
  }
  if (stats.output) {
    sc.stat <- na.omit(sc.stat)
    sc.sum  <- list(coefs = sc.sum, stats = sc.stat)
  } else{
    return(sc.sum)
  }
}
