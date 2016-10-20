arg.vec <- c(
  "test/H3K4me3_TDH_other/jointProblems",
  "test/H3K4me3_TDH_other/joint.model.RData")

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 2){
  stop("usage: Rscript train_model_joint.R samples_dir model.RData")
}
jointProblems <- normalizePath(arg.vec[1], mustWork=TRUE)
joint.model.RData <- arg.vec[2]

library(PeakSegJoint)

target.tsv.vec <- Sys.glob(file.path(
  jointProblems, "*", "target.tsv"))
cat("Found", length(target.tsv.vec), "target.tsv files for training.\n")

problems.list <- list()
for(target.tsv.i in seq_along(target.tsv.vec)){
  target.tsv <- target.tsv.vec[[target.tsv.i]]
  target.vec <- scan(target.tsv, quiet=TRUE)
  problem.dir <- dirname(target.tsv)
  segmentations.RData <- file.path(problem.dir, "segmentations.RData")
  load(segmentations.RData)
  problems.list[[problem.dir]] <- list(
    features=segmentations$features,
    target=target.vec)
}

target.mat <- unname(t(sapply(problems.list, "[[", "target")))
set.seed(1)
n.observations <- length(problems.list)
n.folds <- 3
fold.vec <- sample(rep(1:n.folds, l=n.observations))

squared.hinge <- function(x){
  ifelse(x<1,(x-1)^2,0)
}

error.loss.list <- list()
for(validation.fold in unique(fold.vec)){
  ##print(validation.fold)
  is.validation <- fold.vec == validation.fold
  is.train <- !is.validation
  train.problems <- problems.list[is.train]
  fit <- IntervalRegressionProblems(train.problems, max.iterations=1e3)
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

if(FALSE){
  
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(validation.fold ~ .)+
  geom_point(aes(-log(regularization), mean.surrogate.loss,
                 color=set.name),
             data=validation.min.simplest)+
  geom_line(aes(-log(regularization), mean.surrogate.loss,
                group=paste(validation.fold, set.name),
                color=set.name),
            data=set.error.loss)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(validation.fold ~ .)+
  geom_line(aes(-log(regularization), mean.incorrect.targets,
                group=paste(validation.fold, set.name),
                color=set.name),
            data=set.error.loss)

}

mean.reg <- mean(validation.min.simplest$regularization)
joint.model <- IntervalRegressionProblems(
  problems.list,
  initial.regularization=mean.reg,
  factor.regularization=NULL)

save(joint.model, problems.list, file=joint.model.RData)


