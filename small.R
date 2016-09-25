library(coseg)
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
library(data.table)
loss <- fread(paste0("cat ", problem.dir, "/*_loss.tsv"))
setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status"))
n.bases <- with(one.sample, sum(chromEnd-chromStart))
exp.bases <- rep(n.bases, l=nrow(loss))

library(testthat)
test_that("PeakSegFPOP reports correct number of bases", {
  expect_identical(as.integer(loss$bases), as.integer(exp.bases))
})

loss.ord <- loss[order(penalty),]
loss.ord[, log.penalty := log(penalty)]
best <- loss.ord[target.vec[1] < log.penalty & log.penalty < target.vec[2],]

test_that("Target interval contains one model with 1 peak", {
  expect_identical(as.integer(best$peaks), 1L)
})
