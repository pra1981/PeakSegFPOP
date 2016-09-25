bedGraph.base <-
  "H3K36me3_TDH_immune_McGill0001_chr10_51448845_125869472.bedGraph"
if(!file.exists(bedGraph.base)){
  xz.base <- paste0(bedGraph.base, ".xz")
  xz.file <- system.file(
    "data",
    xz.base,
    package="cosegData")
  file.copy(xz.file, getwd())
  system(paste("unxz", xz.base))
}

penalty.str <- "7400.04974500218"
cmd <- paste("./PeakSegFPOP", bedGraph.base, penalty.str)
status <- system(cmd)

library(testthat)
test_that("FPOP computes correct number of segments", {
  segments.bed <- paste0(
    bedGraph.base, "_penalty=", penalty.str, "_segments.bed")
  segs <- read.table(segments.bed)
  expect_equal(nrow(segs), 1159)
})    

