arg.vec <- "test/H3K36me3_AM_immune_McGill0106_chr16_46385801_88389383"
arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 1){
  stop("usage: Rscript compute_features.R data_dir/sample_dir/problems/problem_dir")
}
problem.dir <- arg.vec[1]

library(coseg)

problem.features(problem.dir)
