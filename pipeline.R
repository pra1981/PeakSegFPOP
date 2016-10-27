arg.vec <- "test/input"

arg.vec <- commandArgs(trailingOnly=TRUE)

set.dir <- normalizePath(arg.vec, mustWork=TRUE)

system.or.stop <- function(cmd){
  cat(cmd, "\n")
  code <- system(cmd)
  if(code != 0){
    stop("non-zero exit code ", code)
  }
}

Rscript <- function(...){
  code <- sprintf(...)
  stopifnot(length(code)==1)
  if(grepl("'", code)){
    print(code)
    stop("there can not be any ' in code")
  }
  sprintf("Rscript -e '%s'", code)
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
  target.cmd <- Rscript('coseg::problem.target("%s")', sample.dir)
  system.or.stop(target.cmd)
})

## Train single-sample model.
train.cmd <- paste("Rscript train_model.R", set.dir)
system.or.stop(train.cmd)

## Single-sample prediction and peak clustering, one job for each
## problem.
sh.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems.bed.sh"))
mclapply.or.stop(sh.vec, function(sh){
  predict.cmd <- paste("bash", sh)
  system.or.stop(predict.cmd)
})

## Compute target intervals for multi-sample problems... does not take
## much time, TODO combine with train_model_joint.R step
labels.tsv.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems", "*", "labels.tsv"))
mclapply.or.stop(labels.tsv.vec, function(labels.tsv){
  target.cmd <- Rscript(
    'PeakSegJoint::problem.joint.target("%s")',
    dirname(labels.tsv))
  system.or.stop(target.cmd)
})
## Train joint model.
train.cmd <- paste("Rscript train_model_joint.R", set.dir)
system.or.stop(train.cmd)

## Joint prediction.
joint.dir.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems", "*"))
joint.model.RData <- file.path(set.dir, "joint.model.RData")
mclapply.or.stop(joint.dir.vec, function(joint.dir), {
  predict.cmd <- Rscript(
    'coseg::problem.predict("%s", "%s")',
    joint.model.RData, joint.dir)
  system.or.stop(predict.cmd)
})

## Summarize peak predictions on a web page.
final.cmd <- paste("Rscript plot_all.R", set.dir)
system.or.stop(final.cmd)
