arg.vec <- "test/H3K4me3_TDH_other/jointProblems/chr4:88923952-88935469"
arg.vec <- "test/H3K4me3_TDH_other/jointProblems/chr11:110160402-110172255"
arg.vec <- "test/input/problems/chr10:38868835-39154935/jointProblems/chr10:39124681-39126535"

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 1){
  stop("usage: Rscript compute_coverage_target.R data_dir/sample_dir/problems/problem_dir")
}
jointProblem.dir <- normalizePath(arg.vec[1], mustWork=TRUE)

PeakSegJoint::problem.joint.target(jointProblem.dir)
