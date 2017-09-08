####
#### R code for Ushio et al.
#### "Fluctuating interaction network and time-varying stability of a natural fish community"
#### No.S3 Twin surrogate
#### (Thiel et al. 2006 Europhys. Lett. 75:535-541)
#### (We would like to thank Dr. Shin-ichiro Nakayama [National Research Institute of Fisheries Science] for his advice on the twin surrogate)
####


# Convert matrix of thresholds into binary 0-1 values
Binalize <- function(e, threshold) {
  if (e > threshold) {
    return(1)
  } else{
    return(0)
  }
}

DistOneZero <- function(m, method = "norm", s = 0.875) {
  # Calculate the maximum norm
  if (method == "norm") {
    dist.mat <-
      apply(m, 1, function(row2) {
        apply(m, 1, function(row1) {
          max(abs(row1 - row2))
        })
      })
  } else if (method == "euclid") {
    dist.mat <- as.matrix(dist(m))
  }

  # Identify 100*s % quantlie value of dist.mat
  d.threshold <- quantile(dist.mat, s)

  # Replace values lower than d.threshold with 0
  # Replace values higher than d.threshold with 1
  rec <- apply(dist.mat, c(1, 2), function(a) Binalize(a, d.threshold))

  # Return binary matrix
  return(rec)
}


# Return next point
PointNext <- function(x, twins, d) {
  if (x == 0) {
    nex <- 0  # Add 0 if the next point is the end point
  } else{
    cand <- c(0, 0)
    if (!is.null(twins[twins[, 1] == x,])) {
      cand <- rbind(cand, twins[twins[, 1] == x,])
      nex <- cand[floor(runif(1, 2, (nrow(cand) + 1))), 2] + 1
    } else{
      nex <- x + 1
    }
  }

  if (nex > ncol(d)) {
    nex <- 0  # Add 0 if the next point is the end point
  }
  nex
}


TwinSurrogate <- function(original,
                          dim,
                          num.iter,
                          tau = 1,
                          s   = 0.875,
                          surrogate.option = "random",  # or "phase_lock"
                          initial.point    = "same_season",  # or "twins"
                          distance.method  = "norm",
                          point.per.year   = 24,
                          s.update         = "on",
                          n.twin.threshold = 10,
                          output.message = F) {
  # Generate time-lag embedding matrix
  if (dim >  1) {
    original_e <- embedd(original, dim, tau)
  }
  if (dim == 1) {
    original_e <- as.matrix(original)
  }
  if (dim <  1) {
    cat("Warning: Embedding dimension should be >= 1")
  }
  # s candidates
  s.seq <- c(0.875, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.87, 0.86, 0.85, 0.84, 0.83, 0.82, 0.81, 0.80)
  s.i <- 1

  # Calculate binary matrix (= recurrence plot)
  repeat {
    d <- DistOneZero(original_e, method = distance.method, s = s)

    # Search twins
    twins <- c()
    for (i in 1:ncol(d)) {
      for (j in 1:ncol(d)) {
        if (all(d[, i] == d[, j])) {
          if (surrogate.option == "phase_lock") {
            if ((i - j) %% point.per.year == 0) {
              twins <- rbind(twins, c(i, j))
            }
          } else{
            twins <- rbind(twins, c(i, j))
          }
        }
      }
    }

    num.twins <- nrow(twins) - nrow(d)
    if (num.twins >= n.twin.threshold) {
      break
    }
    if (s.update == "off") {
      break
    }
    s.i <- s.i + 1
    if (s.i > length(s.seq)) {
      break
    }
    s = s.seq[s.i]
  }

  if (output.message) {
    prop.black <- sum(d) / (ncol(d) * nrow(d))
    print(c("Proportion of 1:", prop.black), quote = F)
    print(c("Number of twins:", num.twins), quote = F)
  }

  # Generate twin surrogates
  surrogate <- as.list(NULL)
  avoid.infinite.loop <- 0
  repeat {
    # Select the initial point of the surrogate
    if (surrogate.option == "random") {
      # Select random initial points
      surr <- sample(1:(ncol(d) - 1), 1)
    } else if (surrogate.option == "phase_lock") {
      if (initial.point == "same_season") {
        # Select the point of the same season
        surr <- sample(seq(1 + point.per.year, ncol(d) - 1, by = point.per.year), 1)
      } else if (initial.point == "twins") {
        # Select twins of the original initial point as the surrogate initial point
        surr <- sample(twins[twins[, 1] == 1, 2], 1)
      }
    } else{
      cat("Warning: specify the correct option!")
    }

    # Search next point
    for (j in 1:(ncol(d) - 1)) {
      repeat {
        nex <- PointNext(surr[length(surr)], twins, d)
        if (surrogate.option == "phase_lock" &&
            initial.point == "same_season") {
          if (nex != (length(surr) + 1)) {
            break
          }
        } else{
          break
        }
      }
      surr <- c(surr, nex)
    }

    # Save the surrogate if it reach to the length of the origial data
    # Not save if the surrogate is short
    if (surr[length(surr)] != 0) {
      surrogate <- c(surrogate, list(original_e[surr,]))
    }

    # Complete the surrogate generation if length(surrogate) == num.iter
    if (length(surrogate) == num.iter) {
      break
    }

    # Avoid infinite loop
    # End cycles if surrogates cannot be generated during > 30*num.iter trials
    avoid.infinite.loop <- avoid.infinite.loop + 1
    if (avoid.infinite.loop > 30 * num.iter) {
      break
    }
  }

  surrogate.one.col <-
    data.frame(matrix(rep(NaN, num.iter * length(original)), ncol = num.iter))
  if (avoid.infinite.loop <= 30 * num.iter) {
    for (i in 1:num.iter) {
      if (dim >= 2) {
        surrogate.one.col[1:(dim - 1), i] <- surrogate[[i]][1, 1:(dim - 1)]
        surrogate.one.col[dim:length(original), i] <- surrogate[[i]][, dim]
      } else {
        surrogate.one.col[dim:length(original), i] <- surrogate[[i]]
      }
    }
  }
  return(surrogate.one.col)
}
