## input: a directory with labels.bed and coverage.bedGraph.

## output: problems directory with sub-directories for each problem
## with labels.

arg.vec <- c("hg19_problems.bed", "labels/H3K36me3_TDH_immune/McGill0001")
arg.vec <- c("hg19_problems.bed", "labels/H3K36me3_AM_immune_folds2-4/McGill0322/")
arg.vec <- commandArgs(trailingOnly=TRUE)
if(length(arg.vec) != 2){
  stop("usage: Rscript problemDB.R problems.bed data_dir/sample_dir/")
}

library(data.table)
library(PeakError)

problems.bed <- normalizePath(arg.vec[1], mustWork=TRUE)
sample.dir <- normalizePath(arg.vec[2], mustWork=TRUE)

problems <- fread(problems.bed)
setnames(problems, c("chrom", "problemStart", "problemEnd"))
problems[, problemStart1 := problemStart + 1L]
problems[, problem.name := sprintf("%s:%d-%d", chrom, problemStart, problemEnd)]
coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")
labels.bed <- file.path(sample.dir, "labels.bed")
labels <- fread(labels.bed)
setnames(labels, c("chrom", "chromStart", "chromEnd", "annotation"))
just.to.check <- PeakError(Peaks(), labels)

labels[, chromStart1 := chromStart + 1L]
problems.dir <- file.path(sample.dir, "problems")
setkey(labels, chrom, chromStart1, chromEnd)
setkey(problems, chrom, problemStart1, problemEnd)
over.dt <- foverlaps(labels, problems, nomatch=0L)
if(nrow(over.dt) < nrow(labels)){
  warning(
    nrow(labels), " lines in ",
    labels.bed, " but only ",
    nrow(over.dt), " labels occur in ",
    problems.bed)
}
labels.by.problem <- split(data.frame(over.dt), over.dt$problem.name)

makeProblem <- function(problem.i){
  problem <- data.frame(problems)[problem.i,]
  problem.dir <- file.path(problems.dir, problem$problem.name)
  cat(sprintf(
    "%4d / %4d %s\n",
    problem.i, nrow(problems), problem.dir))
  dir.create(problem.dir, showWarnings=FALSE, recursive=TRUE)
  problem.bed <- file.path(problem.dir, "problem.bed")
  prob.text <- with(problem, data.frame(
    chrom,
    chromStart=sprintf("%d", problemStart),
    chromEnd=sprintf("%d", problemEnd)))
  write.table(
    prob.text, problem.bed,
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  if(problem$problem.name %in% names(labels.by.problem)){
    problem.labels <- labels.by.problem[[problem$problem.name]]
    prob.lab.bed <- file.path(problem.dir, "labels.bed")
    lab.text <- with(problem.labels, data.frame(
      chrom,
      chromStart=sprintf("%d", chromStart),
      chromEnd=sprintf("%d", chromEnd),
      annotation))
    write.table(
      lab.text,
      prob.lab.bed,
      quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  }
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  sh.file <- paste0(prob.cov.bedGraph, ".sh")
  script.txt <- paste0("#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00
#PBS -A bws-221-ae
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -o ", prob.cov.bedGraph, ".out
#PBS -e ", prob.cov.bedGraph, ".err
#PBS -V                                        
#PBS -N ", problem$problem.name, "
", "Rscript ", normalizePath("compute_coverage_target.R"), " ", problem.dir, "
")
  writeLines(script.txt, sh.file)
}

library(parallel)
nothing <- mclapply(1:nrow(problems), makeProblem)
