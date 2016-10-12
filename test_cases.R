source("test_functions.R")

## Manually create a data set with two chunks from our benchmark.
db.prefix <- "http://cbio.mines-paristech.fr/~thocking/chip-seq-chunk-db/"
set.name <- "H3K4me3_TDH_other"
set.dir <- file.path("test", set.name)
samples.dir <- file.path(set.dir, "samples")
chunk.list <- list(train=c(22, 7), test=24)
chunk.vec <- unlist(chunk.list)
for(chunk.id in chunk.vec){
  chunk.dir <- file.path(set.dir, "chunks", chunk.id)
  chunk.name <- paste0(set.name, "/", chunk.id)
  dir.create(chunk.dir, showWarnings=FALSE, recursive=TRUE)
  counts.RData <- file.path(chunk.dir, "counts.RData")
  if(!file.exists(counts.RData)){
    u <- paste0(db.prefix, chunk.name, "/counts.RData")
    download.file(u, counts.RData)
  }
  load(counts.RData)
  counts.by.sample <- split(counts, counts$sample.id)
  regions.RData <- file.path(chunk.dir, "regions.RData")
  if(!file.exists(regions.RData)){
    u <- paste0(db.prefix, chunk.name, "/regions.RData")
    download.file(u, regions.RData)
  }
  load(regions.RData)
  regions.by.sample <- split(regions, regions$sample.id)
  for(sample.id in names(regions.by.sample)){
    sample.counts <- counts.by.sample[[sample.id]]
    cell.type <- paste(sample.counts$cell.type[1])
    problem.dir <- file.path(
      samples.dir, cell.type, sample.id, "problems", chunk.id)
    dir.create(problem.dir, showWarnings=FALSE, recursive=TRUE)
    sample.regions <- regions.by.sample[[sample.id]]
    sample.counts$chrom <- sample.regions$chrom[1]
    write.table(
      sample.counts[, c("chrom", "chromStart", "chromEnd", "coverage")],
      file.path(problem.dir, "coverage.bedGraph"),
      quote=FALSE,
      sep="\t",
      row.names=FALSE,
      col.names=FALSE)
    problem <- with(sample.counts, data.frame(
      chrom=chrom[1],
      chromStart=chromStart[1],
      chromEnd=max(chromEnd)))
    write.table(
      problem[, c("chrom", "chromStart", "chromEnd")],
      file.path(problem.dir, "problem.bed"),
      quote=FALSE,
      sep="\t",
      row.names=FALSE,
      col.names=FALSE)
    if(chunk.id %in% chunk.list$train){
      write.table(
        sample.regions[, c("chrom", "chromStart", "chromEnd", "annotation")],
        file.path(problem.dir, "labels.bed"),
        quote=FALSE,
        sep="\t",
        row.names=FALSE,
        col.names=FALSE)
      cmd <- paste("Rscript compute_coverage_target.R", problem.dir)
      status <- system(cmd)
      if(status != 0){
        stop("non-zero exit status")
      }
    }
  }
}

loss <- fread("cat test/H3K4me3_TDH_other/samples/kidney/McGill0023/problems/22/*_loss.tsv")
setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status", "mean.intervals", "max.intervals"))
test_that("un-necessary models are not computed", {
  expect_false(31 %in% loss$peaks)
})

model.RData <- file.path(set.dir, "model.RData")
train.cmd <- paste("Rscript train_model.R", samples.dir, model.RData)
system(train.cmd)
test_that("models file created", {
  obj.name.vec <- load(model.RData)
  expect_true("model" %in% obj.name.vec)
})

test.glob <- file.path(samples.dir, "*", "*", "problems", 24)
test.dir.vec <- Sys.glob(test.glob)
segments.glob <- file.path(test.glob, "*_segments.bed")
test.segments.vec <- Sys.glob(segments.glob)
loss.glob <- file.path(test.glob, "*_loss.tsv")
test.loss.vec <- Sys.glob(loss.glob)
peaks.bed.vec <- file.path(test.dir.vec, "peaks.bed")
unlink(peaks.bed.vec)
unlink(test.segments.vec)
unlink(test.loss.vec)
pred.peaks.list <- list()
for(test.dir in test.dir.vec){
  predict.cmd <- paste("Rscript predict_problem.R", model.RData, test.dir)
  system(predict.cmd)
  peaks.bed <- file.path(test.dir, "peaks.bed")
  tryCatch({
    sample.peaks <- fread(peaks.bed)
    setnames(sample.peaks, c("chrom", "chromStart", "chromEnd", "status", "mean"))
    sample.id <- basename(dirname(dirname(test.dir)))
    pred.peaks.list[[test.dir]] <- data.table(sample.id, sample.peaks)
  }, error=function(e){
    ## do nothing
  })
}
test_that("peaks.bed files created", {
  expect_true(all(file.exists(peaks.bed.vec)))
})

## pred.peaks <- do.call(rbind, pred.peaks.list)
## ggplot()+
##   theme_bw()+
##   theme(panel.margin=grid::unit(0, "lines"))+
##   facet_grid(sample.id ~ ., scales="free")+
##   geom_step(aes(chromStart/1e3, coverage),
##             data=counts,
##             color="grey50")+
##   geom_segment(aes(chromStart/1e3, 0,
##                    xend=chromEnd/1e3, yend=0),
##                data=pred.peaks,
##                color="deepskyblue",
##                size=2)

## How to test if peaks are feasible?

## loss <- fread(paste0("cat ", test.glob, "/*_loss.tsv"))
## setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status", "mean.intervals", "max.intervals"))
## test_that("predicted peaks are feasible", {
##   expect_true(all(loss$status=="feasible"))
## })

## this problem has already computed models from 0 to 4 peaks, so we
## should not have to re-run PeakSegFPOP.
problem.dir <- "test/H3K4me3_TDH_other/samples/kidney/McGill0023/problems/7"
loss.files.before <- Sys.glob(file.path("*_loss.tsv"))
predict.cmd <- paste("Rscript predict_problem.R", model.RData, problem.dir)
system(predict.cmd)
loss.files.after <- Sys.glob(file.path("*_loss.tsv"))
test_that("PeakSegFPOP is not run when we already have the solution", {
  expect_equal(length(loss.files.before), length(loss.files.after))
})

## this problem predicts outside the target interval, so we should
## return the closest model inside.
problem.dir <-
  "test/H3K4me3_TDH_other/samples/leukemiaCD19CD10BCells/McGill0267/problems/7"
predict.cmd <- paste("Rscript predict_problem.R", model.RData, problem.dir)
system(predict.cmd)
peaks <- fread(file.path(problem.dir, "peaks.bed"))
test_that("predict model with 2 peaks", {
  expect_equal(nrow(peaks), 2)
})

## Predict for an entire sample using predict_sample.R
sample.dir <- file.path(samples.dir, "skeletalMuscleCtrl", "McGill0036")
sample.pred.cmd <- paste("Rscript predict_sample.R", model.RData, sample.dir)
system(sample.pred.cmd)
test_that("sampleID/peaks.bed file created using predict_sample.R", {
  peaks.bed <- file.path(sample.dir, "peaks.bed")
  peaks <- fread(peaks.bed)
  setnames(peaks, c("chrom", "chromStart", "chromEnd", "status", "mean"))
})

## Predict on some data where we have labels, so we can train the
## joint algo.
labels.bed.vec <- Sys.glob(file.path(
  samples.dir, "*", "*", "problems", "*", "labels.bed"))
for(labels.bed in labels.bed.vec){
  problem.dir <- dirname(labels.bed)
  predict.cmd <- paste("Rscript predict_problem.R", model.RData, problem.dir)
  system(predict.cmd)
}

## Create the scripts that will be used to train the joint algo.
for(chunk.id in chunk.vec){
  cmd <- paste("Rscript create_problems_joint.R", samples.dir, chunk.id)
  system(cmd)
}
jointProblems <- file.path(set.dir, "jointProblems")
target.sh.vec <- Sys.glob(file.path(
  jointProblems, "*", "target.tsv.sh"))
peaks.sh.vec <- Sys.glob(file.path(
  jointProblems, "*", "peaks.bed.sh"))
test_that("more peaks.bed.sh prediction scripts than target.sh training", {
  expect_true(length(target.sh.vec) < length(peaks.sh.vec))
})
test_that("there are some labeled problems", {
  expect_true(0 < length(target.sh.vec))
})

## Compute target intervals.
for(target.sh in target.sh.vec){
  target.cmd <- paste("bash", target.sh)
  system(target.cmd)
}
target.tsv.vec <- sub("[.]sh$", "", target.sh.vec)
test_that("target intervals computed", {
  expect_true(all(file.exists(target.tsv.vec)))
})

## train joint model.
joint.model.RData <- file.path(set.dir, "joint.model.RData")
train.joint.cmd <- paste(
  "Rscript train_model_joint.R", jointProblems, joint.model.RData)
system(train.joint.cmd)
test_that("joint.model.RData created", {
  expect_true(file.exists(joint.model.RData))
})

## Joint prediction.
for(peaks.sh in peaks.sh.vec){
  predict.cmd <- paste("bash", peaks.sh)
  system(predict.cmd)
}
peaks.bed.vec <- sub("[.]sh$", "", peaks.sh.vec)
test_that("joint peaks.bed files created", {
  expect_true(all(file.exists(peaks.bed.vec)))
})

## Longer test for target interval search.
data(H3K36me3_AM_immune_McGill0002_chunk1, package="cosegData")
writeProblem(H3K36me3_AM_immune_McGill0002_chunk1, problem.dir)
test.cmd <- paste("Rscript compute_features.R", problem.dir)
system(test.cmd)
f.row <- read.table(
  file.path(problem.dir, "features.tsv"),
  header=TRUE,
  check.names=FALSE)
test_that("features are computed", {
  expect_equal(nrow(f.row), 1)
})
test.cmd <- paste("Rscript compute_coverage_target.R", problem.dir)
system(test.cmd)
cat.cmd <- paste0("cat ", problem.dir, "/*loss.tsv")
loss <- fread(cat.cmd)
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

target.tsv <- file.path(problem.dir, "target.tsv")
target.vec <- scan(target.tsv, quiet=TRUE)
## 9.18313518325901 12.4640792600544

test_that("target interval is around min error", {
  best <- loss[target.vec[1] < log.penalty & log.penalty < target.vec[2],]
  max.peaks <- max(best$peaks)
  expect_identical(as.integer(max.peaks), 95L)
  min.peaks <- min(best$peaks)
  expect_identical(as.integer(min.peaks), 12L)
})

data("H3K4me3_XJ_immune_chunk1", package="coseg")
sample.id <- "McGill0106"
H3K4me3_XJ_immune_chunk1$count <- H3K4me3_XJ_immune_chunk1$coverage
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
one.sample <- by.sample[[sample.id]]
one.sample$chrom <- "chr1"
n.bases <- with(one.sample, sum(chromEnd-chromStart))
labels.list <- list(
  peakStart=data.frame(
    chrom="chr1",
    chromStart=20007000,
    chromEnd=20009000,
    annotation="peakStart"),
  peaks=data.frame(
    chrom="chr1",
    chromStart=20007000,
    chromEnd=20009000,
    annotation="peaks"),
  noPeaks=data.frame(
    chrom="chr1",
    chromStart=20005000,
    chromEnd=20007000,
    annotation="noPeaks"))

results.list <- list()
samples.dir <- file.path("test", "H3K4me3_TDH_other", "samples")
for(labels.name in names(labels.list)){
  sample.dir <- file.path(samples.dir, "notype", labels.name)
  labels <- labels.list[[labels.name]]
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
  sh.path <- file.path(problem.dir, "target.tsv.sh")
  bash.cmd <- paste("bash", sh.path)
  system(bash.cmd)
  target.tsv <- file.path(problem.dir, "target.tsv")
  target.vec <- scan(target.tsv, quiet=TRUE)
  models.tsv <- file.path(problem.dir, "target_models.tsv")
  models <- fread(models.tsv)
  loss <- fread(paste0("cat ", problem.dir, "/*_loss.tsv"))
  setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status", "mean.intervals", "max.intervals"))
  test_that("sum of loss rows equals target_models rows", {
    expect_equal(nrow(loss), nrow(models))
  })
  test_that("PeakSegFPOP reports correct number of bases", {
    exp.bases <- rep(n.bases, l=nrow(loss))
    expect_identical(as.integer(loss$bases), as.integer(exp.bases))
  })
  loss.ord <- loss[order(penalty),]
  loss.ord[, log.penalty := log(penalty)]
  results.list[[labels.name]] <-
    loss.ord[target.vec[1] <= log.penalty & log.penalty <= target.vec[2],]
}

test_that("Target interval contains one model with 1 peak for peakStart label", {
  expect_equal(results.list$peakStart$peaks, 1)
})

test_that("Target interval contains one model with 1 peak for peaks label", {
  expect_equal(results.list$peaks$peaks, 1)
})

test_that("biggest min error model for noPeaks label has 1 peak", {
  expect_equal(max(results.list$noPeaks$peaks), 1)
})

test_that("smallest min error model for noPeaks label has 0 peaks", {
  expect_equal(min(results.list$noPeaks$peaks), 0)
})

## Predict for an entire sample using sampleID/peaks.bed.sh
peaks.bed <- file.path(sample.dir, "peaks.bed")
peaks.bed.sh <- paste0(peaks.bed, ".sh")
sample.pred.cmd <- paste("bash", peaks.bed.sh)
system(sample.pred.cmd)
test_that("sampleID/peaks.bed file created using peaks.bed.sh", {
  peaks <- fread(peaks.bed)
  setnames(peaks, c("chrom", "chromStart", "chromEnd", "status", "mean"))
})

## Overlapping labels test.
labels <- data.frame(
  chrom="chr1",
  chromStart=c(7, 8),
  chromEnd=c(9, 10),
  annotation="peaks")
sample.dir <- file.path("test", "overlapping")
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

test_that("overlapping labels is an error", {
  status <- system(test.cmd)
  expect_false(status == 0)
})

## Test prediction.
peaks.bed <- file.path(problem.dir, "peaks.bed")
peaks.sh.path <- paste0(peaks.bed, ".sh")
peaks.cmd <- paste("bash", peaks.sh.path)
system(peaks.cmd)
peaks <- fread(peaks.bed)
test_that("peaks.bed has 1 peak predicted", {
  expect_equal(nrow(peaks), 1)
})

