arg.vec <- "test/demo"

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 1){
  stop("usage: Rscript train_model.R data_dir")
}
data.dir <- normalizePath(arg.vec[1], mustWork=TRUE)
samples.dir <- file.path(data.dir, "samples")
model.RData <- file.path(data.dir, "model.RData")

library(data.table)
library(ggplot2)

## For performing K-fold cross-validation in parallel.
if(require(future)){
  plan(multiprocess)
}

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
model <- if(nrow(features) < 10){
  IntervalRegressionUnregularized(
    features[, c("log.quartile.100%", "log.data")], targets)
}else{
  IntervalRegressionCV(
    features, targets, verbose=0,
    initial.regularization=1e-4,
    min.observations=nrow(features),
    reg.type=ifelse(nrow(features) < 20, "1sd", "min(mean)"))
}

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
correct.peaks <- correct.targets[!grepl("Input", problem.dir), {
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

## The size.model is constructed by fitting a normal distribution to
## the log10(bases) values for all the peaks in the correctly
## predicted target intervals in the training data. PeakSegFPOP
## typically gives some peaks which are much larger or smaller than
## the mean, and these may cause problems for the peak clustering
## step. So we use the size.model to exclude peaks which are much
## larger or smaller than the mean. The times parameter is the number
## of standard deviations past log10(mean) which are allowed. Larger
## values of times mean that more peaks will be included in the
## separate peak prediction step (times=Inf means that no peaks will
## be excluded). Smaller values of times mean that fewer peaks will be
## included in the separate peak prediction step (which increases the
## risk of false negatives).
times <- 1
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

## Plot the size model and limits.
log10.bases.grid <- correct.peaks[, seq(
  min(log10.bases), max(log10.bases), l=100)]
normal.dens <- data.table(
  log10.bases=log10.bases.grid,
  prob=size.model[, dnorm(log10.bases.grid, mean, sd)])
base.labels <- size.model[, {
  log10.bases <- c(lower.lim, mean, upper.lim)
  data.table(
    log10.bases,
    hjust=c(1, 0.5, 0),
    label=scales::comma(round(10^log10.bases)))
}]
size.plot <- ggplot()+
  geom_histogram(
    aes(
      log10.bases, ..density..),
    data=correct.peaks)+
  geom_vline(
    aes(xintercept=mean),
    data=size.model,
    size=1, color="red")+
  PeakSegJoint::geom_tallrect(
    aes(xmin=lower.lim, xmax=upper.lim),
    data=size.model, fill="red")+
  geom_line(
    aes(log10.bases, prob),
    data=normal.dens, color="red", size=1)+
  geom_text(
    aes(log10.bases, 0, label=label, hjust=hjust),
    data=base.labels, vjust=1)

cat("Writing model to", model.RData, "\n")
save(
  model, features, targets,
  size.model, size.plot, correct.peaks,
  file=model.RData)

