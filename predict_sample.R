arg.vec <- c(
  "test/H3K4me3_TDH_other/model.RData",
  "test/H3K4me3_TDH_other/samples/McGill0036")

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 2){
  stop("usage: Rscript predict_problem.R model.RData sample_dir")
}
model.RData <- normalizePath(arg.vec[1], mustWork=TRUE)
sample.dir <- normalizePath(arg.vec[2], mustWork=TRUE)

library(data.table)
library(coseg)

problem.bed.vec <- Sys.glob(file.path(sample.dir, "problems", "*", "problem.bed"))
peaks.list <- list()
for(problem.i in seq_along(problem.bed.vec)){
  problem.bed <- problem.bed.vec[[problem.i]]
  cat(sprintf("%4d / %4d %s\n", problem.i, length(problem.bed.vec), problem.bed))
  problem.dir <- dirname(problem.bed)
  peaks.list[[problem.bed]] <- problem.predict(problem.dir, model.RData)
}
peaks <- do.call(rbind, peaks.list)

peaks.bed <- file.path(sample.dir, "peaks.bed")
write.table(
  peaks,
  peaks.bed,
  quote=FALSE,
  col.names=FALSE,
  row.names=FALSE,
  sep="\t")
