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
  expect_equal(nrow(segs), 1157)
})

data(H3K36me3_AM_immune_McGill0002_chunk1, package="cosegData")

problem.dir <- "testProblem"
dir.create(problem.dir, showWarnings=FALSE, recursive=TRUE)
data.list <- H3K36me3_AM_immune_McGill0002_chunk1
file.list <- list(
  coverage=file.path(problem.dir, "coverage.bedGraph"),
  labels=file.path(problem.dir, "labels.bed"),
  problem=file.path(problem.dir, "problem.bed"))
data.list$problem <- with(data.list$coverage, data.frame(
  chrom=chrom[1],
  chromStart=min(chromStart),
  chromEnd=max(chromEnd)))
for(name in names(file.list)){
  df <- data.list[[name]]
  df$chromStart <- sprintf("%d", df$chromStart)
  df$chromEnd <- sprintf("%d", df$chromEnd)
  path <- file.list[[name]]
  write.table(df, path, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

test.cmd <- paste("Rscript compute_coverage_target.R", problem.dir)
system(test.cmd)

test_that("target interval is around min error", {
  
})
