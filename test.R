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

library(data.table)
cat.cmd <- paste0("cat ", problem.dir, "/*loss.tsv")
loss <- fread(cat.cmd)
target.tsv <- file.path(problem.dir, "target.tsv")
target.vec <- scan(target.tsv, quiet=TRUE)
setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status"))
loss[, log.penalty := log(penalty)]

##         penalty  peaks fp fn
##  1:      0.0000 137445  3  0
##  2:    186.2364  13435  3  0
##  3:   1786.3988   1205  2  1
##  4:   3307.1989    419  1  1
##  5:   5484.3420    188  1  1
##  6:   8154.7107    115  1  1
##  7:   9171.8038    101  1  1
##  8:   9433.6006     97  1  1
##  9:   9632.1352     96  1  1
## 10:   9857.0286     95  0  1
## 11:  10502.1280     89  0  1
## 12:  14058.0582     65  0  1
## 13: 202611.5901     14  0  1
## 14: 253616.7421     12  0  1
## 15: 261430.8798     11  0  2
## 16: 274116.0809      9  0  2
## 17: 368444.9776      6  0  2
## 18: 758054.0189      3  0  2
## 19:         Inf      0  0  2
## tdhock@recycled:~/PeakSegFPOP(master*)$

## 9.18313518325901 12.4640792600544
test_that("target interval is around min error", {
  best <- loss[target.vec[1] < log.penalty & log.penalty < target.vec[2],]
  max.peaks <- max(best$peaks)
  expect_identical(as.integer(max.peaks), 95L)
  min.peaks <- min(best$peaks)
  expect_identical(as.integer(min.peaks), 12L)
})
