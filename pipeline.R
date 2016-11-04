arg.vec <- "test/demo"

arg.vec <- commandArgs(trailingOnly=TRUE)

set.dir <- normalizePath(arg.vec, mustWork=TRUE)

system.or.stop <- function(cmd){
  cat(cmd, "\n")
  code <- system(cmd)
  if(code != 0){
    stop("non-zero exit code ", code)
  }
}

## First convert labels.
convert.cmd <- paste("Rscript convert_labels.R", set.dir)
system.or.stop(convert.cmd)

## Create problems for each sample.
create.cmd <- paste("Rscript create_problems_all.R", set.dir)
system.or.stop(create.cmd)

library(coseg)

## Compute target interval for each problem.
samples.dir <- file.path(set.dir, "samples")
labels.bed.vec <- Sys.glob(file.path(
  samples.dir, "*", "*", "problems", "*", "labels.bed"))
mclapply.or.stop(labels.bed.vec, function(labels.bed){
  sample.dir <- dirname(labels.bed)
  coseg::problem.target(sample.dir)
})

## Train single-sample model.
train.cmd <- paste("Rscript train_model.R", set.dir)
system.or.stop(train.cmd)

## Single-sample prediction and peak clustering, one job for each
## problem.
problem.dir.vec <- Sys.glob(file.path(set.dir, "problems", "*"))
for(problem.dir in problem.dir.vec){
  problem.name <- basename(problem.dir)
  create.cmd <- paste(
    "bash",
    file.path(problem.dir, "jointProblems.bed.sh"))
  ## This includes mclapply across samples.
  system.or.stop(create.cmd)
}

## Compute target intervals for multi-sample problems, then learn a
## penalty function for joint peak prediction.
PeakSegJoint::problem.joint.targets.train(set.dir)

## Joint prediction.
joint.dir.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems", "*"))
mclapply.or.stop(joint.dir.vec, PeakSegJoint::problem.joint.predict)

## Summarize peak predictions on a web page.
peaks.tsv.sh <- file.path(set.dir, "peaks_matrix.tsv.sh")
final.cmd <- paste("bash", peaks.tsv.sh)
system.or.stop(final.cmd)
