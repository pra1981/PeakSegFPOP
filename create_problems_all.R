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

arg.vec <- "test/input"
arg.vec <- commandArgs(trailingOnly=TRUE)
if(length(arg.vec) != 1){
  stop("usage: Rscript create_problems_all.R data_dir")
}

library(data.table)
library(PeakError)

data.dir <- normalizePath(arg.vec[1], mustWork=TRUE)
problems.bed <- file.path(data.dir, "problems.bed")
samples.dir <- file.path(data.dir, "samples")
model.RData <- file.path(data.dir, "model.RData")

problems <- fread(problems.bed)
setnames(problems, c("chrom", "problemStart", "problemEnd"))
problems[, problemStart1 := problemStart + 1L]
problems[, problem.name := sprintf(
  "%s:%d-%d", chrom, problemStart, problemEnd)]
cat(
  "Read ", nrow(problems),
  " problems from ", problems.bed,
  "\n", sep="")

sample.dir.vec <- Sys.glob(file.path(samples.dir, "*", "*"))
for(sample.i in seq_along(sample.dir.vec)){
  sample.dir <- sample.dir.vec[[sample.i]]
  problems.dir <- file.path(sample.dir, "problems")
  coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")
  labels.bed <- file.path(sample.dir, "labels.bed")
  labels <- tryCatch({
    labels <- fread(labels.bed)
    setnames(labels, c("chrom", "chromStart", "chromEnd", "annotation"))
    labels
  }, error=function(e){
    cat("No labels in", labels.bed, "\n")
    data.table()
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
    script.txt <- paste0(PBS.header, "
#PBS -o ", target.tsv, ".out
#PBS -e ", target.tsv, ".err
#PBS -N Target", problem$problem.name, "
", "Rscript ", normalizePath("compute_coverage_target.R", mustWork=TRUE), " ", problem.dir, "
")
    writeLines(script.txt, sh.file)
    ## Script for peaks.
    peaks.bed <- file.path(problem.dir, "peaks.bed")
    sh.file <- paste0(peaks.bed, ".sh")
    script.txt <- paste0(PBS.header, "
#PBS -o ", peaks.bed, ".out
#PBS -e ", peaks.bed, ".err
#PBS -N Predict", problem$problem.name, "
", "Rscript ", normalizePath("predict_problem.R", mustWork=TRUE), " ",
model.RData, " ", problem.dir, " 
")
    writeLines(script.txt, sh.file)
  }
  cat(sprintf(
    "Writing %4d / %4d samples %d labels %d labeled problems %s\n",
    sample.i, length(sample.dir.vec),
    nrow(labels), length(labels.by.problem),
    problems.dir))
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
}

## Now write data_dir/problems/*/jointProblems.bed.sh
cat(
  "Writing ", nrow(problems),
  " jointProblems.bed.sh scripts.",
  "\n", sep="")
for(problem.i in 1:nrow(problems)){
  problem <- problems[problem.i,]
  prob.dir <- file.path(data.dir, "problems", problem$problem.name)
  dir.create(prob.dir, showWarnings=FALSE, recursive=TRUE)
  jointProblems.bed <- file.path(prob.dir, "jointProblems.bed")
  sh.file <- paste0(jointProblems.bed, ".sh")
  sample.cmd.vec <- paste(
    "Rscript",
    normalizePath("predict_problem.R", mustWork=TRUE), 
    model.RData,
    file.path(sample.dir.vec, "problems", problem$problem.name))
  script.txt <- paste0(PBS.header, "
#PBS -o ", jointProblems.bed, ".out
#PBS -e ", jointProblems.bed, ".err
#PBS -N P", problem$problem.name, "
", paste(sample.cmd.vec, collapse="\n"), " 
Rscript ",
normalizePath("create_problems_joint.R", mustWork=TRUE),
" ", samples.dir, " ", problem$problem.name, "
")
  writeLines(script.txt, sh.file)
}
