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
samples.dir <- file.path(set.dir, "samples")
sample.dir.vec <- Sys.glob(file.path(samples.dir, "*", "*"))
problems.bed <- file.path(set.dir, "problems.bed")
for(sample.dir in sample.dir.vec){
  create.cmd <- paste("Rscript create_problems.R", problems.bed, sample.dir)
  system.or.stop(create.cmd)
}

## Compute target interval for each problem.
labels.bed.vec <- Sys.glob(file.path(
  samples.dir, "*", "*", "problems", "*", "labels.bed"))
for(labels.bed in labels.bed.vec){
  sample.dir <- dirname(labels.bed)
  create.cmd <- paste("Rscript compute_coverage_target.R", sample.dir)
  system.or.stop(create.cmd)
}

## Train single-sample model.
model.RData <- file.path(set.dir, "model.RData")
train.cmd <- paste("Rscript train_model.R", samples.dir, model.RData)
system.or.stop(train.cmd)

## Single-sample prediction, one job for each problem. Or would it be
## better to create one job for each sample? Or a few jobs for each
## sample?
## library(data.table)
## problems <- fread("hg19_problems.bed")
## setnames(problems, c("chrom", "chromStart", "chromEnd"))
## problems[, bases := chromEnd-chromStart]
## problems[, cum.bases := cumsum(bases)]
## n.jobs <- length(labels.bed.vec)
## bases.per.job <- sum(problems$bases)/n.jobs
## problems[, job.id := cum.bases %/% bases.per.job]
## problems[, list(bases=sum(bases)), by=job.id]
## problem.dir.vec <- Sys.glob(file.path(
##   samples.dir, "*", "*", "problems", "*"))
## for(problem.dir in problem.dir.vec){
##   predict.cmd <- paste("Rscript predict_problem.R", model.RData, problem.dir)
##   system.or.stop(predict.cmd)
## }
for(sample.dir in sample.dir.vec){
  predict.cmd <- paste("Rscript predict_sample.R", model.RData, sample.dir)
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

