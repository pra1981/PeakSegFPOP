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

## First convert labels.
convert.cmd <- paste("Rscript convert_labels.R", set.dir)
system.or.stop(convert.cmd)

## Create problems for each sample.
create.cmd <- paste("Rscript create_problems_all.R", set.dir)
system.or.stop(create.cmd)

## Compute target interval for each problem.
samples.dir <- file.path(set.dir, "samples")
labels.bed.vec <- Sys.glob(file.path(
  samples.dir, "*", "*", "problems", "*", "labels.bed"))
for(labels.bed in labels.bed.vec){
  sample.dir <- dirname(labels.bed)
  target.cmd <- paste("Rscript compute_coverage_target.R", sample.dir)
  system.or.stop(target.cmd)
}

## Train single-sample model.
model.RData <- file.path(set.dir, "model.RData")
train.cmd <- paste("Rscript train_model.R", set.dir)
system.or.stop(train.cmd)

## Single-sample prediction, one job for each problem.
sh.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems.bed.sh"))
for(sh in sh.vec){
  predict.cmd <- paste("bash", sh)
  system.or.stop(predict.cmd)
}

## Create problems for multi-sample segmentation (only depends on peak
## predictions in single-sample problems for which we have labels).
problem.dir.vec <- grep(sample.dir, dirname(labels.bed.vec), value=TRUE)
for(problem.dir in problem.dir.vec){
  problem.name <- basename(problem.dir)
  create.joint.cmd <- paste(
    "Rscript create_problems_joint.R", samples.dir, problem.name)
  system.or.stop(create.joint.cmd)
}

