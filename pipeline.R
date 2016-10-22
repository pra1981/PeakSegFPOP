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
train.cmd <- paste("Rscript train_model.R", set.dir)
system.or.stop(train.cmd)

## Single-sample prediction and peak clustering, one job for each
## problem.
sh.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems.bed.sh"))
for(sh in sh.vec){
  predict.cmd <- paste("bash", sh)
  system.or.stop(predict.cmd)
}

## Compute target intervals for multi-sample problems... how to
## parallelize?
labels.tsv.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems", "*", "labels.tsv"))
for(labels.tsv in labels.tsv.vec){
  target.cmd <- paste("Rscript compute_joint_target.R", dirname(labels.tsv))
  system.or.stop(target.cmd)
}

## Train joint model.
train.cmd <- paste("Rscript train_model_joint.R", set.dir)
system.or.stop(train.cmd)

## Joint prediction.
joint.dir.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems", "*"))
joint.model.RData <- file.path(set.dir, "joint.model.RData")
for(joint.dir in joint.dir.vec){
  predict.cmd <- paste(
    "Rscript predict_problem_joint.R",
    joint.model.RData, joint.dir)
  system.or.stop(predict.cmd)
}

## Plots.
chunk.dir.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "chunks", "*"))
for(chunk.dir in chunk.dir.vec){
  plot.cmd <- paste(
    "Rscript plot_chunk.R",
    chunk.dir)
  system.or.stop(plot.cmd)
}

## TODO: make this a separate script, summarize other info like
## predicted peaks, specific peaks, etc.
figure.png.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "chunks", "*", "figure-predictions-zoomout.png"))
relative.vec <- sub("/", "", sub(set.dir, "", figure.png.vec))
a.vec <- sprintf('
<a href="%s">
  <img src="%s" />
</a>
', relative.vec, sub("-zoomout", "", relative.vec))
writeLines(a.vec, file.path(set.dir, "index.html"))
