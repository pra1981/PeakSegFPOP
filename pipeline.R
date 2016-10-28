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
sample.dir.vec <- Sys.glob(file.path(samples.dir, "*", "*"))
model.RData <- file.path(set.dir, "model.RData")
for(problem.dir in problem.dir.vec){
  problem.name <- basename(problem.dir)
  mclapply.or.stop(sample.dir.vec, function(sample.dir){
    coseg::problem.predict(
      file.path(sample.dir, "problems", problem.name),
      model.RData)
  })
  create.cmd <- paste("Rscript create_problems_joint.R", samples.dir, problem.name)
  system.or.stop(create.cmd)
}

## Compute target intervals for multi-sample problems... does not take
## much time, TODO combine with train_model_joint.R step
labels.tsv.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems", "*", "labels.tsv"))
mclapply.or.stop(labels.tsv.vec, function(labels.tsv){
  PeakSegJoint::problem.joint.target(dirname(labels.tsv))
})
## Train joint model.
train.cmd <- paste("Rscript train_model_joint.R", set.dir)
system.or.stop(train.cmd)

## Joint prediction.
joint.dir.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems", "*"))
joint.model.RData <- file.path(set.dir, "joint.model.RData")
mclapply.or.stop(joint.dir.vec, function(joint.dir){
  PeakSegJoint::problem.joint.predict(joint.model.RData, joint.dir)
})

## Summarize peak predictions on a web page.
final.cmd <- paste("Rscript plot_all.R", set.dir)
system.or.stop(final.cmd)
