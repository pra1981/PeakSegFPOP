arg.vec <- c(
  "test/H3K4me3_TDH_other")

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 1){
  stop("usage: Rscript train_model_joint.R data_dir")
}
data.dir <- normalizePath(arg.vec[1], mustWork=TRUE)

PeakSegJoint::problem.joint.train(data.dir)

