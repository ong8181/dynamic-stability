####
#### R code for Ushio et al.
#### "Fluctuating interaction network and time-varying stability of a natural fish community"
#### Configuration
####

config <- list()

config$kBestE.Method   <- "MAE"
config$kBestE.Range    <- 1:24
config$kBestE.RangeStr <-
  paste0("E", config$kBestE.Range[1], "to",
         config$kBestE.Range[length(config$kBestE.Range)])
config$kFishNameFileEncoding <- "Shift-JIS"
config$kMacFont <- "HiraMaruProN-W4"
config$kPdfFont <- "Japan1GothicBBB"
config$kMaxCore <- 0
config$kRndSeed <- 2430

if (!is.null(config$kRndSeed)) {
  set.seed(config$kRndSeed)
}
