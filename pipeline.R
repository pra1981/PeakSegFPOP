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
sample.dir.vec <- Sys.glob(file.path(set.dir, "samples", "*", "*"))
problems.bed <- file.path(set.dir, "problems.bed")
for(sample.dir in sample.dir.vec){
  create.cmd <- paste("Rscript create_problems.R", problems.bed, sample.dir)
  system.or.stop(create.cmd)
}

## Compute target interval for each problem.
labels.bed.vec <- Sys.glob(file.path(set.dir, "samples", "*", "*", "problems", "*", "labels.bed"))
for(labels.bed in labels.bed.vec){
  sample.dir <- dirname(labels.bed)
  create.cmd <- paste("Rscript compute_coverage_target.R", sample.dir)
  system.or.stop(create.cmd)
}
