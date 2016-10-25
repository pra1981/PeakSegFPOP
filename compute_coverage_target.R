## inputs: problem directory with problem.bed.

## outputs: if the problem directory does not have coverage.bedGraph,
## it will be created via intersectBed. If labels.bed is present, we
## also create target.bed.
arg.vec <- "test/PeakSegJoint-two/samples/tcell/McGill0322/problems/chr11:96437584-134946516"
arg.vec <- "test/H3K36me3_AM_immune_McGill0002_chunk1"
arg.vec <- "labels/H3K36me3_AM_immune_folds2-4/McGill0002/problems/chr1:3995268-13052998"
arg.vec <- "labels/H3K36me3_TDH_immune/McGill0001/problems/chr11:96437584-134946516"
arg.vec <- "test/noPeaks/problems/chr1:17175658-29878082"
arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 1){
  stop("usage: Rscript compute_coverage_target.R data_dir/sample_dir/problems/problem_dir")
}
problem.dir <- normalizePath(arg.vec[1], mustWork=TRUE)

library(coseg)
library(data.table)
library(PeakError)

t.info <- problem.target(problem.dir)

