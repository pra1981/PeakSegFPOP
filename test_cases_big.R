source("test_functions.R")

problem.dir <- "test/H3K36me3_AM_immune_McGill0107_chr13_19020000_86760324.RData"
obj.name <- data(H3K36me3_AM_immune_McGill0107_chr13_19020000_86760324, package="cosegData")
data.list <- get(obj.name)
data.list$coverage <- data.list$coverage[1:3100000,]
writeProblem(data.list, problem.dir)
coverage.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
test.cmd <- paste("./PeakSegFPOP", coverage.bedGraph, data.list$penalty)
system(test.cmd)

problem.dir <- "test/H3K36me3_AM_immune_McGill0029_chr16_46385801_88389383"
obj.name <- data(H3K36me3_AM_immune_McGill0029_chr16_46385801_88389383, package="cosegData")
data.list <- get(obj.name)
data.list$coverage <- data.list$coverage[1:1000000,]
writeProblem(data.list, problem.dir)
coverage.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
test.cmd <- paste("./PeakSegFPOP", coverage.bedGraph, data.list$penalty)
system(test.cmd)

problem.dir <- "test/H3K36me3_AM_immune_McGill0106_chr16_46385801_88389383"
data(H3K36me3_AM_immune_McGill0106_chr16_46385801_88389383, package="cosegData")
writeProblem(H3K36me3_AM_immune_McGill0106_chr16_46385801_88389383, problem.dir)
test.cmd <- paste("Rscript compute_coverage_target.R", problem.dir)
system(test.cmd)

cat.cmd <- paste0("cat ", problem.dir, "/*loss.tsv")
loss <- fread(cat.cmd)
setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status", "mean.intervals", "max.intervals"))
test_that("un-helpful models are not computed", {
  expect_false(any(1:6 %in% loss$peaks))
})
target.tsv <- file.path(problem.dir, "target.tsv")
target.vec <- scan(target.tsv, quiet=TRUE)
test_that("target interval computed", {
  expect_equal(length(target.vec), 2)
})
loss[, log.penalty := log(penalty)]
optimal <- loss[target.vec[1] < log.penalty & log.penalty < target.vec[2], ]
test_that("target interval is correct", {
  expect_equal(range(optimal$peaks), c(43, 610))
})

problem.dir <- "test/H3K36me3_AM_immune_McGill0079_chr3_60000_66170270"
data(H3K36me3_AM_immune_McGill0079_chr3_60000_66170270, package="cosegData")
writeProblem(H3K36me3_AM_immune_McGill0079_chr3_60000_66170270, problem.dir)
test.cmd <- paste("Rscript compute_coverage_target.R", problem.dir)
system(test.cmd)

cat.cmd <- paste0("cat ", problem.dir, "/*loss.tsv")
loss <- fread(cat.cmd)
setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status", "mean.intervals", "max.intervals"))
test_that("un-helpful models are not computed", {
  expect_false(60 %in% loss$peaks)
  expect_false(2084 %in% loss$peaks)
})

target.tsv <- file.path(problem.dir, "target.tsv")
test_that("target interval computed", {
  target.vec <- scan(target.tsv, quiet=TRUE)
  expect_equal(length(target.vec), 2)
})

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
cmd <- paste("PeakSegFPOP", bedGraph.base, penalty.str)
status <- system(cmd)

test_that("FPOP computes correct number of segments", {
  segments.bed <- paste0(
    bedGraph.base, "_penalty=", penalty.str, "_segments.bed")
  segs <- read.table(segments.bed)
  expect_equal(nrow(segs), 1157)
})

