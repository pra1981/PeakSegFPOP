## input: a directory with labels.bed and coverage.bedGraph.

## output: problems directory with sub-directories for each problem
## with labels.

arg.vec <- c("hg19_problems.bed", "labels/H3K36me3_TDH_immune/McGill0001")
arg.vec <- commandArgs(trailingOnly=TRUE)
if(length(arg.vec) != 2){
  stop("usage: Rscript problemDB.R problems.bed data_dir/sample_dir/")
}

library(data.table)
problems.bed <- arg.vec[1]
sample.dir <- arg.vec[2]

problems <- fread(problems.bed)
setnames(problems, c("chrom", "problemStart", "problemEnd"))
problems[, problemStart1 := problemStart + 1L]
problems[, problem.name := paste0(chrom, ":", problemStart, "-", problemEnd)]
coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")
labels.bed <- file.path(sample.dir, "labels.bed")
labels <- fread(labels.bed)
setnames(labels, c("chrom", "chromStart", "chromEnd", "annotation"))
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
  pname <- names(labels.by.problem)[[problem.i]]
  problem.dir <- file.path(problems.dir, pname)
  problem.labels <- labels.by.problem[[problem.i]]
  problem <- problem.labels[1, c("chrom", "problemStart", "problemEnd")]
  dir.create(problem.dir, showWarnings=FALSE, recursive=TRUE)
  problem.bed <- file.path(problem.dir, "problem.bed")
  write.table(
    problem, problem.bed,
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  prob.lab.bed <- file.path(problem.dir, "labels.bed")
  write.table(
    problem.labels[, c("chrom", "chromStart", "chromEnd", "annotation")],
    prob.lab.bed,
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  cat(sprintf(
    "%4d / %4d %s\n",
    problem.i, length(labels.by.problem), problem.dir))
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  if(!file.exists(prob.cov.bedGraph)){
    ## Use intersectBed to avoid memory problems.
    cmd <- paste(
      "intersectBed -sorted",
      "-a", coverage.bedGraph,
      "-b", problem.bed,
      ">", prob.cov.bedGraph)
    cat(cmd, "\n")
    system(cmd)
  }
}

library(parallel)
nothing <- mclapply(seq_along(labels.by.problem), makeProblem)
