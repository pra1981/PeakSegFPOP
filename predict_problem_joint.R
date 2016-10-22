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

library(PeakSegJoint)
library(data.table)

load(joint.model.RData)

converted <- problem.joint(jointProblem.dir)
log.penalty <- joint.model$predict(converted$features)
stopifnot(length(log.penalty)==1)
selected <- subset(
  converted$modelSelection,
  min.log.lambda < log.penalty & log.penalty < max.log.lambda)

pred.dt <- if(selected$peaks == 0){
  data.table()
}else{
  pred.df <- subset(converted$peaks, peaks==selected$peaks)
  chrom <- paste(converted$coverage$chrom[1])
  with(pred.df, data.table(
    chrom,
    chromStart,
    chromEnd,
    name=paste0(sample.group, "/", sample.id),
    mean))
}
peaks.bed <- file.path(jointProblem.dir, "peaks.bed")
cat("Writing ",
    nrow(pred.dt), " peaks to ",
    peaks.bed,
    "\n", sep="")
write.table(
  pred.dt, peaks.bed,
  quote=FALSE,
  sep="\t",
  col.names=FALSE,
  row.names=FALSE)
