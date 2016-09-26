library(testthat)
library(data.table)
library(coseg)

data(H3K36me3_AM_immune_McGill0002_chunk1, package="cosegData")

problem.dir <- "testProblem"
unlink(problem.dir, recursive=TRUE)
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

cat.cmd <- paste0("cat ", problem.dir, "/*loss.tsv")
loss <- fread(cat.cmd)
target.tsv <- file.path(problem.dir, "target.tsv")
target.vec <- scan(target.tsv, quiet=TRUE)
setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status", "mean.intervals", "max.intervals"))
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

data("H3K4me3_XJ_immune_chunk1")
sample.id <- "McGill0106"
H3K4me3_XJ_immune_chunk1$count <- H3K4me3_XJ_immune_chunk1$coverage
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
one.sample <- by.sample[[sample.id]]
one.sample$chrom <- "chr1"
labels <- data.frame(
  chrom="chr1",
  chromStart=20007000,
  chromEnd=20009000,
  annotation="peakStart")
sample.dir <- "labels/small/McGill0106"
unlink(sample.dir, recursive=TRUE)
dir.create(sample.dir, showWarnings=FALSE, recursive=TRUE)
coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")
write.table(
  one.sample[, c("chrom", "chromStart", "chromEnd", "coverage")],
  coverage.bedGraph,
  quote=FALSE, sep="\t",
  row.names=FALSE, col.names=FALSE)
labels.bed <- file.path(sample.dir, "labels.bed")
write.table(
  labels, labels.bed,
  quote=FALSE, sep="\t",
  row.names=FALSE, col.names=FALSE)

test.cmd <- paste("Rscript create_problems.R hg19_problems.bed", sample.dir)
system(test.cmd)
labels.bed <- Sys.glob(paste0(sample.dir, "/problems/*/labels.bed"))
problem.dir <- dirname(labels.bed)
sh.path <- file.path(problem.dir, "coverage.bedGraph.sh")
bash.cmd <- paste("bash", sh.path)
system(bash.cmd)
target.tsv <- file.path(problem.dir, "target.tsv")
target.vec <- scan(target.tsv, quiet=TRUE)
loss <- fread(paste0("cat ", problem.dir, "/*_loss.tsv"))
setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status", "mean.intervals", "max.intervals"))
n.bases <- with(one.sample, sum(chromEnd-chromStart))
exp.bases <- rep(n.bases, l=nrow(loss))

test_that("PeakSegFPOP reports correct number of bases", {
  expect_identical(as.integer(loss$bases), as.integer(exp.bases))
})

loss.ord <- loss[order(penalty),]
loss.ord[, log.penalty := log(penalty)]
best <- loss.ord[target.vec[1] < log.penalty & log.penalty < target.vec[2],]

test_that("Target interval contains one model with 1 peak", {
  expect_identical(as.integer(best$peaks), 1L)
})
