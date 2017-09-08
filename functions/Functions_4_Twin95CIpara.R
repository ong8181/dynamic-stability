####
#### R code for Ushio et al.
#### "Fluctuating interaction network and time-varying stability of a natural fish community"
#### No.S4 Parallel version of function to calculate 95%CI using twin surrogate
####

# Function to calculate 95% CI for the surrogate (Parallel version)
ParSurCI <- function(target.ts,
                     cause.ts,
                     surrogate.ts,
                     lib.parms = c(11, 15, 24, 40, 80, 120, 160, 200, 240, length(target.ts)),
                     surrogate = "effect",
                     E.range   = 1:24,
                     tp        = 0) {
  E.tar <- BestErEDM(target.ts, E = E.range)
  E.cau <- BestErEDM(cause.ts, E = E.range)
  lib.m <- length(target.ts)
  lib.size.s <- lib.parms
  surrogate.sum <- data.frame(L = lib.size.s)

  # Do CCM for the surrogate data
  surrogate.all <-
    pforeach(
      i      = 1:ncol(surrogate.ts),
      .c     = cbind,
      .cores = config$kMaxCore,
      .seed  = config$kRndSeed
    )({
      if (surrogate == "effect") {
        target.sur <- surrogate.ts[, i]
        block      <- cbind(target.sur, cause.ts)
      } else if (surrogate == "cause") {
        cause.sur <-  surrogate.ts[, i]
        block     <- cbind(target.ts, cause.sur)
      }
      m <- nrow(block)

      ccm.tar.cau <-
        ccm(
          block,
          E         = E.tar,
          lib_sizes = lib.size.s,
          silent    = T,
          tp        = tp
        )
      ccm.m <- ccm_means(ccm.tar.cau)[, c('lib_size', 'rho')]

      rhos.tmp <- ccm.m$rho
    })

  if (!is.null(surrogate.all)) {
    surrogate.sum <- as.data.frame(cbind(surrogate.sum, surrogate.all[1:length(lib.parms), 1:ncol(surrogate.ts)]))
  }

  # calculate 95% CI
  ccm.sur.ci <-
    data.frame(
      L = lib.size.s,
      lower95 = NA,
      upper95 = NA,
      lower99 = NA,
      upper99 = NA
    )
  csa <- surrogate.sum
  n.s <- length(surrogate.ts)

  for(j in 1:nrow(surrogate.sum)){
    upper95 <- quantile(as.numeric(csa[j,2:(n.s+1)]), 0.975)
    lower95 <- quantile(as.numeric(csa[j,2:(n.s+1)]), 0.025)
    ccm.sur.ci[j,'upper95'] <- upper95
    ccm.sur.ci[j,'lower95'] <- lower95
    upper99 <- quantile(as.numeric(csa[j,2:(n.s+1)]), 0.995)
    lower99 <- quantile(as.numeric(csa[j,2:(n.s+1)]), 0.005)
    ccm.sur.ci[j,'upper99'] <- upper99
    ccm.sur.ci[j,'lower99'] <- lower99
  }

  result <- list(rho = ccm.sur.ci)
  return(result)
}
