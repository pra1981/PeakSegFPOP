arg.vec <- c(
  "test/H3K4me3_TDH_other/joint.model.RData",
  "test/H3K4me3_TDH_other/jointProblems/chr4:88923952-88935469")
arg.vec <- c(
  "test/input/joint.model.RData",
  "test/input/problems/chr10:18024675-38818835/jointProblems/chr10:33693003-33974944")

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 2){
  stop("usage: Rscript predict_problem_joint.R model.RData data_dir/jointProblems/problem_dir")
}
joint.model.RData <- normalizePath(arg.vec[1], mustWork=TRUE)
jointProblem.dir <- normalizePath(arg.vec[2], mustWork=TRUE)

PeakSegJoint::problem.joint.predict(joint.model.RData, jointProblem.dir)

