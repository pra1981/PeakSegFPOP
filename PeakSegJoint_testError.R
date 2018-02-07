arg.vec <- "labels/H3K36me3_TDH_immune"

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec)!=1){
  stop("usage: Rscript testError.R path/to/project_dir")
}
data.dir <- normalizePath(arg.vec, mustWork=TRUE)

library(penaltyLearning)
library(PeakSegJoint)

joint.model.RData <- file.path(data.dir, "joint.model.RData")
target.tsv.vec <- Sys.glob(file.path(
  data.dir, "problems", "*", "jointProblems", "*", "target.tsv"))
cat("Found", length(target.tsv.vec), "target.tsv files for training.\n")
problems.list <- list()
target.mat.list <- list()
feature.mat.list <- list()
for(target.tsv.i in seq_along(target.tsv.vec)){
  target.tsv <- target.tsv.vec[[target.tsv.i]]
  target.vec <- scan(target.tsv, quiet=TRUE)
  problem.dir <- dirname(target.tsv)
  cat(sprintf("%4d / %4d %s\n", target.tsv.i, length(target.tsv.vec), problem.dir))
  segmentations.RData <- file.path(problem.dir, "segmentations.RData")
  load(segmentations.RData)
  if(any(is.finite(target.vec))){
    problems.list[[problem.dir]] <- list(
      features=segmentations$features,
      target=target.vec)
    target.mat.list[[problem.dir]] <- target.vec
    feature.mat.list[[problem.dir]] <- colSums(segmentations$features)
  }
}
feature.mat <- do.call(rbind, feature.mat.list)
target.mat <- do.call(rbind, target.mat.list)

n.folds <- 4
set.seed(1)
outer.fold.vec <- sample(rep(1:n.folds, l=nrow(target.mat)))

CVProb <- function(problems.list){
  n.observations <- length(problems.list)
  n.folds <- ifelse(n.observations < 10, 3, 5)
  fold.vec <- sample(rep(1:n.folds, l=n.observations))
  squared.hinge <- function(x){
    ifelse(x<1,(x-1)^2,0)
  }
  error.loss.list <- list()
  for(validation.fold in unique(fold.vec)){
    ##cat(sprintf("%d"))
    is.validation <- fold.vec == validation.fold
    is.train <- !is.validation
    train.problems <- problems.list[is.train]
    fit <- IntervalRegressionProblems(
      train.problems, max.iterations=1e4, verbose=0)
    for(problem.i in seq_along(problems.list)){
      problem <- problems.list[[problem.i]]
      log.penalty.vec <- fit$predict(problem$features)
      too.lo <- log.penalty.vec < problem$target[1]
      too.hi <- problem$target[2] < log.penalty.vec
      left.term <- squared.hinge(log.penalty.vec-problem$target[1])
      right.term <- squared.hinge(problem$target[2]-log.penalty.vec)
      error.loss.list[[paste(validation.fold, problem.i)]] <- data.table(
        validation.fold,
        surrogate.loss=left.term + right.term,
        incorrect.targets=too.lo | too.hi,
        set.name=ifelse(is.train[problem.i], "train", "validation"),
        regularization=fit$regularization.vec,
        problem.i)
    }
  }
  error.loss <- do.call(rbind, error.loss.list)
  set.error.loss <- error.loss[, list(
    mean.surrogate.loss=mean(surrogate.loss),
    mean.incorrect.targets=mean(incorrect.targets)
    ), by=.(validation.fold, set.name, regularization)]
  validation.min <- set.error.loss[set.name=="validation", {
    .SD[mean.surrogate.loss==min(mean.surrogate.loss),]
  }, by=validation.fold]
  validation.min.simplest <- validation.min[, {
    .SD[regularization==max(regularization),]
  }, by=validation.fold]
  mean.reg <- mean(validation.min.simplest$regularization)
  IntervalRegressionProblems(
    problems.list,
    initial.regularization=mean.reg,
    factor.regularization=NULL,
    verbose=0)
}

test.error.list <- list()
for(test.fold in 1:n.folds){
  is.test <- outer.fold.vec == test.fold
  is.train <- !is.test
  train.problems <- problems.list[is.train]
  train.features <- feature.mat[is.train,]
  train.targets <- target.mat[is.train,]
  fit.mat <- IntervalRegressionCV(train.features, train.targets)
  fit.prob <- CVProb(train.problems)
  pred.list <- list(
    mat=fit.mat$predict(feature.mat[is.test,]),
    prob=sapply(which(is.test), function(i)fit.prob$predict(problems.list[[i]]$features)))
  test.targets <- target.mat[is.test,]
  n.test <- sum(is.test)
  for(model.name in names(pred.list)){
    pred.vec <- pred.list[[model.name]]
    is.hi <- test.targets[,2] < pred.vec 
    is.lo <- pred.vec < test.targets[,1]
    incorrect.targets <- sum(is.hi|is.lo)
    test.error.list[[paste(test.fold, model.name)]] <- data.table(
      test.fold, model.name, percent.incorrect=100*incorrect.targets/n.test)
  }
}
test.error <- do.call(rbind, test.error.list)

dcast(test.error, test.fold ~ model.name)
