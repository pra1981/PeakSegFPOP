source("test_functions.R")
writeProblem <- function(...){}

## Test for target interval after infeasible models.
obj.name <- "H3K4me3_TDH_immune_tcell_McGill0007_chr10_60000_17974675"
problem.dir <- file.path("test", obj.name)
data(list=obj.name, package="cosegData")
data.list <- get(obj.name)
writeProblem(data.list, problem.dir)
test.cmd <- Rscript('coseg::problem.target("%s")', problem.dir)
system(test.cmd)
target.vec <- scan(file.path(problem.dir, "target.tsv"), quiet=TRUE)
models <- fread(file.path(problem.dir, "target_models.tsv"))
models.above <- models[target.vec[2] < log(penalty), ]
n.infeasible <- sum(models.above$status == "infeasible")
test_that("only feasible models above upper penalty limit", {
  expect_equal(n.infeasible, 0)
})

## Test for not computing a useless model that does not help find the target interval.
obj.name <- "H3K4me3_TDH_immune_McGill0005_chr1_17175658_29878082"
problem.dir <- file.path("test", obj.name)
data(list=obj.name, package="cosegData")
data.list <- get(obj.name)
writeProblem(data.list, problem.dir)
test.cmd <- Rscript('coseg::problem.target("%s")', problem.dir)
system(test.cmd)
models <- fread(file.path(problem.dir, "target_models.tsv"))
test_that("did not compute infeasible model", {
  expect_false(506 %in% models$peaks)
})

obj.name <- "H3K4me3_TDH_immune_McGill0322_chr3_93504854_194041961"
problem.dir <- file.path("test", obj.name)
data(list=obj.name, package="cosegData")
data.list <- get(obj.name)
writeProblem(data.list, problem.dir)
test.cmd <- Rscript('coseg::problem.target("%s")', problem.dir)
system(test.cmd)
models <- fread(file.path(problem.dir, "target_models.tsv"))
test_that("did not compute feasible limit models", {
  expect_true(all(! (343:344) %in% models$peaks))
})

problem.dir <- "test/H3K36me3_AM_immune_McGill0027_chr4_75452279_191044276"
obj.name <- data(H3K36me3_AM_immune_McGill0027_chr4_75452279_191044276, package="cosegData")
data.list <- get(obj.name)
writeProblem(data.list, problem.dir)
test.cmd <- Rscript('coseg::problem.target("%s")', problem.dir)
system(test.cmd)
peaks.with.fp <- 15:22
target.vec <- scan(file.path(problem.dir, "target.tsv"), quiet=TRUE)
loss <- fread(paste0("cat ", problem.dir, "/*_loss.tsv"))
setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status", "mean.intervals", "max.intervals"))
loss[, log.penalty := log(penalty)]
in.target <- loss[target.vec[1] < log.penalty & log.penalty < target.vec[2],]
test_that("target interval contains no fp models", {
  expect_true(all(!peaks.with.fp %in% in.target))
})
test_that("target interval contains no infeasible models", {
  expect_true(all(in.target$status == "feasible"))
})
min.peaks <- min(in.target$peaks)
smaller.present <- (min.peaks-1) %in% loss$peaks
loss.ord <- loss[order(penalty),]
min.peaks.i <- which(loss.ord$peaks==min.peaks)
min.and.next <- loss.ord[c(min.peaks.i, min.peaks.i+1),]
peaks.count.tab <- table(loss$peaks)
min.and.next.counts <- peaks.count.tab[paste(min.and.next$peaks)]
two.present <- any(1 < min.and.next.counts)
test_that("upper limit computed precisely", {
  expect_true(two.present || smaller.present)
})

problem.dir <- "test/H3K36me3_AM_immune_McGill0106_chr16_46385801_88389383"
data(H3K36me3_AM_immune_McGill0106_chr16_46385801_88389383, package="cosegData")
writeProblem(H3K36me3_AM_immune_McGill0106_chr16_46385801_88389383, problem.dir)
test.cmd <- Rscript('coseg::problem.target("%s")', problem.dir)
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
  expect_equal(range(optimal$peaks), c(43, 113))
})

problem.dir <- "test/H3K36me3_AM_immune_McGill0079_chr3_60000_66170270"
data(H3K36me3_AM_immune_McGill0079_chr3_60000_66170270, package="cosegData")
writeProblem(H3K36me3_AM_immune_McGill0079_chr3_60000_66170270, problem.dir)
test.cmd <- Rscript('coseg::problem.target("%s")', problem.dir)
system(test.cmd)
cat.cmd <- paste0("cat ", problem.dir, "/*loss.tsv")
loss <- fread(cat.cmd)
setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status", "mean.intervals", "max.intervals"))
test_that("un-helpful models are not computed", {
  expect_false(2084 %in% loss$peaks)
  expect_false(618 %in% loss$peaks)
})
target.tsv <- file.path(problem.dir, "target.tsv")
test_that("target interval computed", {
  target.vec <- scan(target.tsv, quiet=TRUE)
  expect_equal(length(target.vec), 2)
})

## Below we have some data for which the minimum cost was not computed
## correctly (problems in C++ code not R scripts).

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

problem.dir <- "test/H3K36me3_AM_immune_McGill0107_chr13_19020000_86760324"
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

