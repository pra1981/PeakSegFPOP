## Edit the following definition to reflect your cluster
## configuration.
PBS.header <- "#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -A bws-221-ae
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -V"

## input: a directory with labels.bed and coverage.bedGraph.

## output: problems directory with sub-directories for each problem
## with labels.

arg.vec <- c("hg19_problems.bed", "test/H3K4me3_TDH_other/samples/peakStart")
arg.vec <- c("hg19_problems.bed", "labels/H3K36me3_AM_immune_folds2-4/McGill0322/")
arg.vec <- commandArgs(trailingOnly=TRUE)
if(length(arg.vec) != 2){
  stop("usage: Rscript problemDB.R problems.bed data_dir/sample_dir/")
}

library(data.table)
library(PeakError)

Rscript <- function(...){
  code <- sprintf(...)
  stopifnot(length(code)==1)
  if(grepl("'", code)){
    print(code)
    stop("there can not be any ' in code")
  }
  sprintf("Rscript -e '%s'", code)
}

problems.bed <- normalizePath(arg.vec[1], mustWork=TRUE)
sample.dir <- normalizePath(arg.vec[2], mustWork=TRUE)
group.dir <- dirname(sample.dir)
samples.dir <- dirname(group.dir)
data.dir <- dirname(samples.dir)
model.RData <- file.path(data.dir, "model.RData")
problems.dir <- file.path(sample.dir, "problems")

problems <- fread(problems.bed)
setnames(problems, c("chrom", "problemStart", "problemEnd"))
problems[, problemStart1 := problemStart + 1L]
problems[, problem.name := sprintf(
  "%s:%d-%d", chrom, problemStart, problemEnd)]
coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")

labels.bed <- file.path(sample.dir, "labels.bed")
labels <- tryCatch({
  labels <- fread(labels.bed)
  setnames(labels, c("chrom", "chromStart", "chromEnd", "annotation"))
  labels
}, error=function(e){
  cat("No labels in", labels.bed, "\n")
  list()
})

labels.by.problem <- if(length(labels)){
  just.to.check <- PeakError(Peaks(), labels)
  labels[, chromStart1 := chromStart + 1L]
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
  split(data.frame(over.dt), over.dt$problem.name)
}

makeProblem <- function(problem.i){
  problem <- data.frame(problems)[problem.i,]
  problem.dir <- file.path(problems.dir, problem$problem.name)
  ## cat(sprintf(
  ##   "%4d / %4d %s\n",
  ##   problem.i, nrow(problems), problem.dir))
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
  ## Script for coverage.
  target.tsv <- file.path(problem.dir, "target.tsv")
  sh.file <- paste0(target.tsv, ".sh")
  target.code <- sprintf('coseg::problem.target("%s")', problem.dir)
  target.cmd <- sprintf("Rscript -e '%s'", target.code)
  script.txt <- paste0(PBS.header, "
#PBS -o ", target.tsv, ".out
#PBS -e ", target.tsv, ".err
#PBS -N Target", problem$problem.name, "
", target.cmd, "
")
  writeLines(script.txt, sh.file)
  ## Script for peaks.
  peaks.bed <- file.path(problem.dir, "peaks.bed")
  sh.file <- paste0(peaks.bed, ".sh")
  predict.cmd <- Rscript(
    'coseg::problem.predict("%s", "%s")',
    model.RData, problem.dir)
  script.txt <- paste0(PBS.header, "
#PBS -o ", peaks.bed, ".out
#PBS -e ", peaks.bed, ".err
#PBS -N Predict", problem$problem.name, "
", predict.cmd, " 
")
  writeLines(script.txt, sh.file)
}

cat("Creating ", nrow(problems), " problems in ", problems.dir, ".\n", sep="")
nothing <- lapply(1:nrow(problems), makeProblem)

## Script for peaks on the whole sample.
peaks.bed <- file.path(sample.dir, "peaks.bed")
sh.file <- paste0(peaks.bed, ".sh")
sample.id <- basename(sample.dir)
script.txt <- paste0(PBS.header, "
#PBS -o ", peaks.bed, ".out
#PBS -e ", peaks.bed, ".err
#PBS -N Predict", sample.id, "
", "Rscript ", normalizePath("predict_sample.R", mustWork=TRUE), " ",
model.RData, " ", sample.dir, " 
")
writeLines(script.txt, sh.file)
