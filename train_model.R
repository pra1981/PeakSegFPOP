arg.vec <- "test/demo"

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 1){
  stop("usage: Rscript train_model.R data_dir")
}
data.dir <- normalizePath(arg.vec[1], mustWork=TRUE)
samples.dir <- file.path(data.dir, "samples")
model.RData <- file.path(data.dir, "model.RData")

library(data.table)

## For performing K-fold cross-validation in parallel.
library(doParallel)
registerDoParallel()

glob.str <- file.path(
  samples.dir, "*", "*", "problems", "*", "target.tsv")
cat("Searching for", glob.str, "files for training.\n")
target.tsv.vec <- Sys.glob(glob.str)
cat("Found", length(target.tsv.vec), "target.tsv files for training.\n")

features.list <- list()
targets.list <- list()
for(target.tsv.i in seq_along(target.tsv.vec)){
  target.tsv <- target.tsv.vec[[target.tsv.i]]
  problem.dir <- dirname(target.tsv)
  features.tsv <- file.path(problem.dir, "features.tsv")
  if(!file.exists(features.tsv)){
    cat(sprintf("%4d / %4d Computing %s\n", target.tsv.i, length(target.tsv.vec), features.tsv))
    problem.features(problem.dir)
  }
  features.list[[problem.dir]] <- fread(features.tsv)
  targets.list[[problem.dir]] <- scan(target.tsv, quiet=TRUE)
}
features <- as.matrix(do.call(rbind, features.list))
targets <- do.call(rbind, targets.list)

library(penaltyLearning)
set.seed(1)
model <- IntervalRegressionCV(
  features, targets, verbose=0,
  min.observations=nrow(features))

cat("Learned regularization parameter and weights:\n")
print(model$pred.param.mat)
pred.log.penalty <- as.numeric(model$predict(features))
pred.dt <- data.table(
  problem.dir=dirname(target.tsv.vec),
  too.lo=as.logical(pred.log.penalty < targets[,1]),
  lower.limit=targets[,1],
  pred.log.penalty,
  upper.limit=targets[,2],
  too.hi=as.logical(targets[,2] < pred.log.penalty))
pred.dt[, status := ifelse(
  too.lo, "low",
  ifelse(too.hi, "high", "correct"))]
correct.targets <- pred.dt[status=="correct", ]
correct.peaks <- correct.targets[, {
  target_models.tsv <- file.path(problem.dir, "target_models.tsv")
  target.models <- fread(target_models.tsv)
  closest <- target.models[which.min(abs(log(penalty)-pred.log.penalty)),]
  coverage.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  segments.bed <- paste0(
    coverage.bedGraph, "_penalty=", closest$penalty, "_segments.bed")
  segs <- fread(segments.bed)
  setnames(segs, c("chrom", "chromStart", "chromEnd", "status", "mean"))
  segs[status=="peak", ]
}, by=problem.dir]
correct.peaks[, bases := chromEnd-chromStart]
correct.peaks[, log10.bases := log10(bases)]
size.model <- correct.peaks[, list(
  mean=mean(log10.bases),
  sd=sd(log10.bases)
)]
times <- 2
size.model[, upper.lim := mean + times*sd]
size.model[, lower.lim := mean - times*sd]
size.model[, upper.bases := 10^(upper.lim)]
size.model[, lower.bases := 10^(lower.lim)]

cat("Train errors:\n")
pred.dt[, list(targets=.N), by=status]

## To check if we are extrapolating when predictions are made later,
## we save the range of the data 
model$train.feature.ranges <- apply(
  features[, model$pred.feature.names, drop=FALSE], 2, range)

cat("Writing model to", model.RData, "\n")
save(model, features, targets, size.model, correct.peaks, file=model.RData)

