source("test_functions.R")

options(warn=2)
test_that("downloaded hg38 problems is the same as stored", {
  system("Rscript downloadProblems.R hg38 hg38_test.bed")
  diff.lines <- system("diff hg38_test.bed hg38_problems.bed", intern=TRUE)
  expect_equal(length(diff.lines), 0)
})
options(warn=0)

status <- system("PeakSegFPOP notIntBad.bedGraph 0.1")
test_that("PeakSegFPOP fails for non-integer data", {
  expect_false(status == 0)
})

status <- system("PeakSegFPOP missing.bedGraph 0.1")
test_that("PeakSegFPOP fails for missing data", {
  expect_false(status == 0)
})

## PeakSegJoint example data.
exampleData <- system.file("exampleData", package="PeakSegJoint")
bigwig.vec <- Sys.glob(file.path(exampleData, "*", "*.bigwig"))
label.files.list <- list(
  overlapping=c(
    "overlapping_labels.txt", "manually_annotated_region_labels.txt"),
  one=c("manually_annotated_region_labels.txt"),
  two=c("manually_annotated_region_labels.txt", "other_labels.txt"))
for(set.name in names(label.files.list)){
  set.dir <- file.path("test", paste0("PeakSegJoint-", set.name))
  for(bigwig.old in bigwig.vec){
    group.dir.old <- dirname(bigwig.old)
    sample.group <- basename(group.dir.old)
    base.old <- basename(bigwig.old)
    sample.id <- sub("[.]bigwig$", "", base.old)
    bigwig.new <- file.path(
      set.dir, "samples", sample.group, sample.id, "coverage.bigWig")
    sample.dir.new <- dirname(bigwig.new)
    dir.create(sample.dir.new, showWarnings=FALSE, recursive=TRUE)
    file.symlink(bigwig.old, bigwig.new)
  }
  label.files <- label.files.list[[set.name]]
  labels.dir <- file.path(set.dir, "labels")
  dir.create(labels.dir, showWarnings=FALSE)
  for(label.file in label.files){
    label.old <- file.path(exampleData, label.file)
    label.new <- file.path(labels.dir, label.file)
    file.symlink(label.old, label.new)
  }
}

## Data set with one bad overlapping label file that should error.
convert.cmd <- "Rscript convert_labels.R test/PeakSegJoint-overlapping"
status <- system(convert.cmd)
test_that("stop for overlapping labels", {
  expect_equal(status, 1)
})

## Data set with one good label file.
set.dir <- file.path("test", "PeakSegJoint-one")
convert.cmd <- paste("Rscript convert_labels.R", set.dir)
status <- system(convert.cmd)
test_that("converting one labels file succeeds", {
  expect_equal(status, 0)
})
labels.bed.vec <- Sys.glob(file.path(
  set.dir, "samples", "*", "*", "labels.bed"))
sample.dir.vec <- dirname(labels.bed.vec)
group.dir.vec <- dirname(sample.dir.vec)
group.name.vec <- basename(group.dir.vec)
test_that("labels.bed files for tcell and bcell samples", {
  expect_equal(sum(group.name.vec=="bcell"), 2)
  expect_equal(sum(group.name.vec=="tcell"), 2)
})
test_that("no labels.bed files for other samples", {
  expect_equal(sum(group.name.vec=="kidney"), 0)
  expect_equal(sum(group.name.vec=="stem"), 0)
  expect_equal(sum(group.name.vec=="skeletalMuscle"), 0)
})
chunk.limits.RData <- file.path(set.dir, "chunk.limits.RData")
test_that("chunk.limits file is created", {
  expect_true(file.exists(chunk.limits.RData))
})

## Data set with two good label files.
set.dir <- file.path("test", "PeakSegJoint-two")
convert.cmd <- paste("Rscript convert_labels.R", set.dir)
status <- system(convert.cmd)
test_that("converting two labels file succeeds", {
  expect_equal(status, 0)
})
labels.bed.vec <- Sys.glob(file.path(
  set.dir, "samples", "*", "*", "labels.bed"))
sample.dir.vec <- dirname(labels.bed.vec)
group.dir.vec <- dirname(sample.dir.vec)
group.name.vec <- basename(group.dir.vec)
test_that("labels.bed files for tcell and bcell samples", {
  expect_equal(sum(group.name.vec=="bcell"), 2)
  expect_equal(sum(group.name.vec=="tcell"), 2)
})
test_that("labels.bed files for other samples", {
  expect_equal(sum(group.name.vec=="kidney"), 1)
  expect_equal(sum(group.name.vec=="stem"), 1)
  expect_equal(sum(group.name.vec=="skeletalMuscle"), 2)
})
immune.lines <- readLines(file.path(
  set.dir, "labels", "manually_annotated_region_labels.txt"))
n.immune <- sum(immune.lines != "")
other.lines <- readLines(file.path(
  set.dir, "labels", "other_labels.txt"))
n.other <- sum(other.lines != "")
n.labels.vec <- rep(NA_integer_, length(labels.bed.vec))
for(labels.bed.i in seq_along(labels.bed.vec)){
  labels.bed <- labels.bed.vec[[labels.bed.i]]
  label.lines <- readLines(labels.bed)
  n.labels.vec[[labels.bed.i]] <- length(label.lines)
}
n.expected.vec <- ifelse(
  group.name.vec %in% c("tcell", "bcell"), n.immune, n.other)
test_that("expected number of lines in each labels.bed file", {
  expect_equal(n.labels.vec, n.expected.vec)
})
chunk.limits.RData <- file.path(set.dir, "chunk.limits.RData")
test_that("chunk.limits file is created", {
  expect_true(file.exists(chunk.limits.RData))
})

## Create problems.
sample.dir <- dirname(labels.bed)
create.cmd <- paste(
  "Rscript create_problems_sample.R hg19_problems.bed", sample.dir)
system(create.cmd)
labels.bed.vec <- Sys.glob(file.path(
  sample.dir, "problems", "*", "labels.bed"))
test_that("at least one labeled problem", {
  expect_gt(length(labels.bed.vec), 0)
})

## Compute target interval.
problem.dir <- dirname(labels.bed.vec[1])
target.tsv <- file.path(problem.dir, "target.tsv")
unlink(target.tsv)
unlink(Sys.glob(file.path(problem.dir, "*_loss.tsv")))
coseg::problem.target(problem.dir)
test_that("target.tsv file created", {
  expect_true(file.exists(target.tsv))
})

## Manually create a data set with two chunks from our benchmark.
db.prefix <- "http://cbio.mines-paristech.fr/~thocking/chip-seq-chunk-db/"
set.name <- "H3K4me3_TDH_other"
set.dir <- file.path("test", set.name)
samples.dir <- file.path(set.dir, "samples")
chunk.list <- list(train=c(22, 7), test=24)
chunk.vec <- unlist(chunk.list)
counts.by.chunk <- list()
regions.by.chunk <- list()
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
  counts.by.chunk[[paste(chunk.id)]] <- counts
  counts.by.sample <- split(counts, counts$sample.id)
  regions.RData <- file.path(chunk.dir, "regions.RData")
  if(!file.exists(regions.RData)){
    u <- paste0(db.prefix, chunk.name, "/regions.RData")
    download.file(u, regions.RData)
  }
  load(regions.RData)
  regions.by.chunk[[paste(chunk.id)]] <- regions
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
      coseg::problem.target(problem.dir)
    }
  }
}

## Should not compute two 0 peaks models in this data set where there
## is a 1 peak model.
models <- fread("test/H3K4me3_TDH_other/samples/skeletalMuscleCtrl/McGill0036/problems/22/target_models.tsv")
test_that("only one model with 0 peaks", {
  expect_equal(sum(models$peaks==0), 1)
})

## The model with 31 peaks does not help to find the target interval,
## so we should not compute it.
models <- fread("cat test/H3K4me3_TDH_other/samples/kidney/McGill0023/problems/22/target_models.tsv")
test_that("un-necessary models are not computed", {
  expect_false(31 %in% models$peaks)
})

## it should be possible to achieve zero errors in each of these
## simple problems.
models <- fread("test/H3K4me3_TDH_other/samples/kidney/McGill0023/problems/22/target_models.tsv")
test_that("target interval includes no errors", {
  expect_true(0 %in% models[, fp + fn])
})
models <- fread("test/H3K4me3_TDH_other/samples/skeletalMuscleCtrl/McGill0037/problems/7/target_models.tsv")
test_that("target interval includes no errors", {
  expect_true(0 %in% models[, fp + fn])
})
models <- fread("test/H3K4me3_TDH_other/samples/skeletalMuscleMD/McGill0016/problems/7/target_models.tsv")
test_that("target interval includes no errors", {
  expect_true(0 %in% models[, fp + fn])
})

model.RData <- file.path(set.dir, "model.RData")
train.cmd <- paste("Rscript train_model.R", set.dir)
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
  coseg::problem.predict(test.dir)
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
coseg::problem.predict(problem.dir)
loss.files.after <- Sys.glob(file.path("*_loss.tsv"))
test_that("PeakSegFPOP is not run when we already have the solution", {
  expect_equal(length(loss.files.before), length(loss.files.after))
})

## Predict for an entire sample using predict_sample.R
sample.dir <- file.path(samples.dir, "skeletalMuscleCtrl", "McGill0036")
sample.pred.cmd <- paste("Rscript predict_sample.R", sample.dir)
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
  coseg::problem.predict(problem.dir)
}

## Create the scripts that will be used to train the joint algo.
for(chunk.id in chunk.vec){
  prob.dir <- file.path(set.dir, "problems", chunk.id)
  cmd <- paste("Rscript create_problems_joint.R", prob.dir)
  system(cmd)
}
target.sh.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems", "*", "target.tsv.sh"))
peaks.sh.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems", "*", "peaks.bed.sh"))
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
train.joint.cmd <- paste("Rscript train_model_joint.R", set.dir)
system(train.joint.cmd)
joint.model.RData <- file.path(set.dir, "joint.model.RData")
test_that("joint.model.RData created", {
  expect_true(file.exists(joint.model.RData))
})

## Joint prediction.
all.peaks.list <- list()
all.problems.list <- list()
for(peaks.sh in peaks.sh.vec){
  predict.cmd <- paste("bash", peaks.sh)
  peaks.bed <- sub("[.]sh$", "", peaks.sh)
  problem.bed <- sub("peaks.bed$", "problem.bed", peaks.bed)
  problem <- fread(problem.bed)
  setnames(problem, c("chrom", "problemStart", "problemEnd", "chunk.id"))
  all.problems.list[[problem.bed]] <- problem
  unlink(peaks.bed)
  system(predict.cmd)
  test_that("joint peaks.bed created", {
    expect_true(all(file.exists(peaks.bed)))
  })
  tryCatch({
    problem.peaks <- fread(peaks.bed)
    setnames(problem.peaks, c(
      "chrom", "chromStart", "chromEnd", "sample.path", "mean"))
    all.peaks.list[[peaks.sh]] <- data.table(
      chunk.id=problem$chunk.id,
      problem.peaks)
  }, error=function(e){
    NULL
  })
}
all.peaks <- do.call(rbind, all.peaks.list)
all.peaks[, sample.id := sub(".*/", "", sample.path)]
all.peaks[, sample.group := sub("/.*", "", sample.path)]
peaks.by.chunk <- split(all.peaks, all.peaks$chunk.id)
all.problems <- do.call(rbind, all.problems.list)
problems.by.chunk <- split(all.problems, all.problems$chunk.id)

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
for(chunk.id in chunk.vec){
  chunk.peaks <- peaks.by.chunk[[paste(chunk.id)]]
  chunk.regions <- regions.by.chunk[[paste(chunk.id)]]
  chunk.regions$sample.group <- chunk.regions$cell.type
  chunk.errors <- PeakErrorSamples(chunk.peaks, chunk.regions)
  if(interactive()){
    first.peaks.list <- list()
    peaks.bed.vec <- Sys.glob(file.path(
      samples.dir, "*", "*", "problems", chunk.id, "peaks.bed"))
    for(peaks.bed in peaks.bed.vec){
      first.peaks.list[[peaks.bed]] <- tryCatch({
        sample.peaks <- fread(peaks.bed)
        chunk.dir <- dirname(peaks.bed)
        problems.dir <- dirname(chunk.dir)
        sample.dir <- dirname(problems.dir)
        sample.id <- basename(sample.dir)
        group.dir <- dirname(sample.dir)
        sample.group <- basename(group.dir)
        setnames(sample.peaks, c(
          "chrom", "chromStart", "chromEnd", "status", "mean"))
        data.table(sample.id, sample.group, sample.peaks)
      }, error=function(e){
        NULL
      })
    }
    first.peaks <- do.call(rbind, first.peaks.list)
    chunk.problems <- problems.by.chunk[[paste(chunk.id)]]
    chunk.counts <- counts.by.chunk[[paste(chunk.id)]]
    gg <- ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(sample.id ~ ., scales="free")+
      scale_fill_manual(values=ann.colors)+
      geom_tallrect(aes(
        xmin=problemStart/1e3,
        xmax=problemEnd/1e3),
        data=chunk.problems,
        fill="grey",
        color="black")+
    geom_tallrect(aes(
      xmin=chromStart/1e3,
      xmax=chromEnd/1e3,
      fill=annotation),
      color="grey",
      alpha=0.5,
      data=chunk.regions)+
    geom_tallrect(aes(
      xmin=chromStart/1e3,
      xmax=chromEnd/1e3,
      linetype=status),
                  color="black",
                  size=1,
                  fill=NA,
      data=chunk.errors)+
    scale_linetype_manual(
      "error type",
      limits=c("correct", 
               "false negative",
               "false positive"),
      values=c(correct=0,
               "false negative"=3,
               "false positive"=1))+
    geom_step(aes(chromStart/1e3, coverage),
              data=chunk.counts,
              color="grey50")+
    geom_segment(aes(chromStart/1e3, 0,
                     xend=chromEnd/1e3, yend=0),
                 data=chunk.peaks,
                 color="deepskyblue",
                 size=2)+
      geom_segment(aes(chromStart/1e3, 0,
                       xend=chromEnd/1e3, yend=0),
                   data=first.peaks)
    print(gg)
  }
  with(chunk.errors, cat(sprintf(
    "chunk=%4d FP=%4d/%4d FN=%4d/%4d\n",
    chunk.id,
    sum(fp), sum(possible.fp),
    sum(fn), sum(possible.tp))))
}

## Longer test for target interval search.
prob.name <- "H3K36me3_AM_immune_McGill0002_chunk1"
data(list=prob.name, package="cosegData")
data.list <- get(prob.name)
problem.dir <- file.path("test", prob.name)
writeProblem(data.list, problem.dir)
test.cmd <- paste("Rscript compute_features.R", problem.dir)
system(test.cmd)
f.row <- read.table(
  file.path(problem.dir, "features.tsv"),
  header=TRUE,
  check.names=FALSE)
test_that("features are computed", {
  expect_equal(nrow(f.row), 1)
})
coseg::problem.target(problem.dir)
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
  test.cmd <- paste("Rscript create_problems_sample.R hg19_problems.bed", sample.dir)
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

test_that("Target interval contains a 1 peak model for peaks label", {
  expect_true(1 %in% results.list$peaks$peaks)
})

test_that("Target interval contains a 0 peak model for noPeaks label", {
  expect_true(0 %in% results.list$noPeaks$peaks)
})

test_that("smallest min error model for noPeaks label has 0 peaks", {
  expect_equal(min(results.list$noPeaks$peaks), 0)
})

## First line of segments.bed should be the last segment, and its
## chromEnd should be the same as the last line of coverage.bedGraph.
segments.glob <- paste0(problem.dir, "/*_segments.bed")
head.cmd <- paste0("head -1 -q ", segments.glob)
last.segments <- fread(head.cmd)
setnames(last.segments, c("chrom", "chromStart", "chromEnd", "status", "mean"))
tail.cmd <- paste("tail -1", coverage.bedGraph)
last.coverage <- fread(tail.cmd)
setnames(last.coverage, c("chrom", "chromStart", "chromEnd", "count"))
test_that("last segments chromEnd == last coverage chromEnd", {
  end.vec <- rep(last.coverage$chromEnd, nrow(last.segments))
  expect_equal(last.segments$chromEnd, end.vec)
})

## Last line of segments.bed should be the first segment, and its
## chromStart should be the same as the first line of coverage.bedGraph.
tail.cmd <- paste0("tail -n 1 -q ", segments.glob)
first.segments <- fread(tail.cmd)
setnames(first.segments, c("chrom", "chromStart", "chromEnd", "status", "mean"))
head.cmd <- paste("head -1", coverage.bedGraph)
first.coverage <- fread(head.cmd)
setnames(first.coverage, c("chrom", "chromStart", "chromEnd", "count"))
test_that("first segments chromStart == first coverage chromStart", {
  start.vec <- rep(first.coverage$chromStart, nrow(first.segments))
  expect_equal(first.segments$chromStart, start.vec)
})

## Make sure create_problems_sample.R works for samples with no labels.
unlink(file.path(sample.dir, "labels.bed"))
peaks.bed <- file.path(sample.dir, "peaks.bed")
peaks.bed.sh <- paste0(peaks.bed, ".sh")
unlink(peaks.bed.sh)
test.cmd <- paste("Rscript create_problems_sample.R hg19_problems.bed", sample.dir)
status <- system(test.cmd)
test_that("script finishes successfully", {
  expect_equal(status, 0)
})
test_that("sampleID/peaks.bed.sh is created", {
  expect_true(file.exists(peaks.bed.sh))
})

## Predict for an entire sample using sampleID/peaks.bed.sh
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
test.cmd <- paste(
  "Rscript create_problems_sample.R hg19_problems.bed", sample.dir)

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


