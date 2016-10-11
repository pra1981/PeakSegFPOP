arg.vec <- c(
  "test/H3K4me3_TDH_other/jointProblems",
  "test/H3K4me3_TDH_other/joint.model.RData")

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 2){
  stop("usage: Rscript train_model_joint.R samples_dir model.RData")
}
jointProblems <- normalizePath(arg.vec[1], mustWork=TRUE)
model.RData <- arg.vec[2]

library(PeakSegJoint)

target.tsv.vec <- Sys.glob(file.path(
  jointProblems, "*", "target.tsv"))
cat("Found", length(target.tsv.vec), "target.tsv files for training.\n")

features.list <- list()
targets.list <- list()
for(target.tsv.i in seq_along(target.tsv.vec)){
  target.tsv <- target.tsv.vec[[target.tsv.i]]
  problem.dir <- dirname(target.tsv)
  features.tsv <- file.path(problem.dir, "features.tsv")
  if(!file.exists(features.tsv)){
    stop("TODO implement joint feature computation")
  }
  features.list[[problem.dir]] <- fread(features.tsv)
  targets.list[[problem.dir]] <- scan(target.tsv, quiet=TRUE)
}
features <- as.matrix(do.call(rbind, features.list))
targets <- do.call(rbind, targets.list)

model <- IntervalRegressionMatrixCV(features, targets, verbose=1)

save(model, features, targets, file=model.RData)


