arg.vec <- c(#predict new data, target not yet computed.
  "test/H3K4me3_TDH_other/model.RData",
  "test/H3K4me3_TDH_other/samples/McGill0036/problems/24")
arg.vec <- c(#predicts inside, already computed.
  "test/H3K4me3_TDH_other/model.RData",
  "test/H3K4me3_TDH_other/samples/McGill0036/problems/22")
arg.vec <- c(#predicts outside target interval
  "test/H3K4me3_TDH_other/model.RData",
  "test/H3K4me3_TDH_other/samples/McGill0267/problems/7")
arg.vec <- c(#predicts infeasible.
  "test/H3K4me3_TDH_other/model.RData",
  "test/H3K4me3_TDH_other/samples/McGill0012/problems/24")

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 2){
  stop("usage: Rscript predict_problem.R model.RData problem_dir")
}
model.RData <- arg.vec[1]
problem.dir <- arg.vec[2]

library(data.table)
library(coseg)

load(model.RData)

features.tsv <- file.path(problem.dir, "features.tsv")
if(!file.exists(features.tsv)){
  cat(sprintf("Computing %s\n", features.tsv))
  problem.features(problem.dir)
}

features <- fread(features.tsv)
feature.mat <- as.matrix(features)
pred.penalty <- as.numeric(exp(model$predict(feature.mat)))
cat(sprintf(
  "Predicting penalty=%f log(penalty)=%f based on %d features.\n",
  pred.penalty,
  log(pred.penalty),
  length(model$pred.feature.names)
  ))

loss.ord <- tryCatch({
  loss <- suppressWarnings(fread(paste0("cat ", problem.dir, "/*_loss.tsv")))
  setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status", "mean.intervals", "max.intervals"))
  loss[, log.penalty := log(penalty)]
  loss[order(-penalty),]
}, error=function(e){
  data.table()
})

## TODO: If we have already computed the target interval and the
## prediction is outside, then we should choose the minimal error
## model which is closest to the predicted penalty.
target.vec <- tryCatch({
  suppressWarnings(scan(file.path(problem.dir, "target.tsv"), quiet=TRUE))
}, error=function(e){
  NULL
})
if(nrow(loss.ord) && length(target.vec)==2){
  cat(sprintf("Target interval %f < log(penalty) < %f\n", target.vec[1], target.vec[2]))
  if(log(pred.penalty) < target.vec[1]){
    pred.penalty <- loss.ord[target.vec[1] < log.penalty, penalty[.N]]
    cat(sprintf(
      "Closest inside target penalty=%f log(penalty)=%f\n",
      pred.penalty, log(pred.penalty)))
  }
  if(target.vec[2] < log(pred.penalty)){
    pred.penalty <- loss.ord[log.penalty < target.vec[2], penalty[1]]
    cat(sprintf(
      "Closest inside target penalty=%f log(penalty)=%f\n",
      pred.penalty, log(pred.penalty)))
  }
}  

## This will be NULL until we find or compute a model that can be used
## for predicted peaks.
pen.str <- NULL

## If two neighboring penalties have already been computed, then we do
## not have to re-run PeakSegFPOP.
if(is.null(pen.str) && 2 <= nrow(loss.ord)){
  is.after <- loss.ord[, penalty < pred.penalty]
  first.after <- which(is.after)[1]
  last.before <- first.after - 1
  smaller.peaks <- loss.ord[last.before, peaks]
  bigger.peaks <- loss.ord[first.after, peaks]
  if(smaller.peaks + 1 == bigger.peaks){
    loss.unique <- unique(loss.ord[, .(peaks, total.cost)])
    exact <- loss.unique[, exactModelSelection(total.cost, peaks, peaks)]
    selected <- subset(exact, min.lambda < pred.penalty & pred.penalty < max.lambda)
    same.peaks <- loss.ord[peaks==selected$peaks, ]
    pen.num <- same.peaks$penalty[1]
    cat(
      "Based on previous computations, penalty of ",
      pred.penalty, " and ", 
      pen.num, " both recover ",
      selected$peaks, " peak",
      ifelse(selected$peaks==1, "", "s"), ".\n",
      sep="")
    pen.str <- paste(pen.num)
  }
}

## TODO: If we have not already computed the target interval, then we
## can run PeakSegFPOP at the predicted penalty value. If the
## resulting model is feasible then we are done. Otherwise, we need to
## compute the target interval to find the biggest feasible model,
## which we return.
if(is.null(pen.str)){
  pen.str <- paste(pred.penalty)
  result <- problem.PeakSegFPOP(problem.dir, pen.str)
  if(result$loss$status=="infeasible"){
    t.info <- problem.target(problem.dir)
    biggest.feasible <- t.info$models[which.min(errors),]
    pen.str <- paste(biggest.feasible$penalty)
  }
}

## compute peaks.
prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
pre <- paste0(prob.cov.bedGraph, "_penalty=", pen.str)
penalty_segments.bed <- paste0(pre, "_segments.bed")
penalty.segs <- fread(penalty_segments.bed)
setnames(penalty.segs, c("chrom","chromStart", "chromEnd", "status", "mean"))
peaks <- penalty.segs[status=="peak", ]
peaks.bed <- file.path(problem.dir, "peaks.bed")
cat(
  "Writing ", peaks.bed,
  " with ", nrow(peaks),
  " peak", ifelse(nrow(peaks)==1, "", "s"),
  " based on ", penalty_segments.bed,
  ".\n", sep="")

write.table(
  peaks,
  peaks.bed,
  quote=FALSE,
  sep="\t",
  col.names=FALSE,
  row.names=FALSE)
