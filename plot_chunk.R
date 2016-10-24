arg.vec <- "test/input/problems/chr10:18024675-38818835/chunks/chr10:33061897-33974942"
arg.vec <- "test/input/problems/chr10:38868835-39154935/chunks/chr10:39098319-39140858"

arg.vec <- commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(data.table)

if(length(arg.vec) != 1){
  stop("usage: Rscript plot_chunk.R project/problems/problemID/chunks/chunkID")
}

chunk.dir <- normalizePath(arg.vec, mustWork=TRUE)
PeakSegJoint::problem.joint.plot(chunk.dir)
